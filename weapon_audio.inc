// Offline Weapon Lab voices. Frequencies are Hz, envelope times are seconds.
#include <stdatomic.h>
enum { LAB_RATE = 44100 };
static bool labPrototype;
static Sound labBeam, labContact, labMotors[32];
static void lab_release_sustain(Sound *sound)
{
    if (sound->frameCount) {
        UnloadSoundAlias(*sound);
        *sound = (Sound){0};
    }
}
static unsigned labTake[LAB_VOICE_COUNT];
static AudioStream labChargeStream;
static atomic_int labChargeLevel;
static float labPhaseA, labPhaseB, labChargeGain;

static float lab_noise(uint32_t *seed)
{
    *seed = *seed * 1664525u + 1013904223u;
    return (*seed >> 8) / 8388608.f - 1;
}

// Bake a diffuse field of reflections into the mono voice so it follows the shot's pan.
static void lab_lightning_reverb(short *pcm, unsigned frames, int take)
{
    enum { EXCITATION = LAB_RATE / 5 };
    float damped[EXCITATION], filtered = 0;
    for (unsigned i = 0; i < EXCITATION; ++i) {
        filtered += .36f * (pcm[i] - filtered);
        damped[i] = filtered;
    }
    float *wet = calloc(frames, sizeof(*wet));
    if (!wet) return;
    uint32_t seed = 713u + take;
    for (unsigned delay = LAB_RATE / 42; delay < frames; delay += 51 + (seed % 97)) {
        float gain = .20f * lab_noise(&seed) * expf(-6.8f * delay / LAB_RATE);
        unsigned count = frames-delay < EXCITATION ? frames-delay : EXCITATION;
        for (unsigned i = 0; i < count; ++i) wet[delay+i] += gain*damped[i];
    }
    for (unsigned i = 0; i < frames; ++i) {
        float fade = fminf(1, (frames-1-i)/(LAB_RATE*.025f));
        pcm[i] = (short)fmaxf(-32767, fminf(32767, pcm[i] + fade*wet[i]));
    }
    free(wet);
}

static Wave lab_voice_wave(int voice, int take)
{
    const float lengths[] = {.40f, .85f, .19f, .21f, .24f, .76f, .28f, .14f};
    unsigned frames = (unsigned)(lengths[voice] * LAB_RATE);
    short *pcm = MemAlloc(frames * sizeof(*pcm));
    if (!pcm) return (Wave){0};
    float low = 0, variation = (take - 1.5f) * .03f;
    float railRumble = 0;
    uint32_t seed = 32171u + voice;
    for (unsigned i = 0; i < frames; i++) {
        float t = (float)i / LAB_RATE;
        if (i == LAB_RATE / 100) seed += take * 7193u;
        float noise = lab_noise(&seed);
        low += .22f * (noise - low);
        float y;
        if (voice == LAB_SHOTGUN) {
            float clack = fmaxf(0, t - (.089f + variation * .16f));
            y = .67f * expf(-23*t) * sinf(2*PI*(58*t + 1.95f*(1-expf(-48*t))))
              + .85f * low * expf(-37*t)
              + .25f * (low + sinf(2*PI*(710+variation*900)*clack)) * expf(-88*clack) * fminf(1,clack*1800);
        } else if (voice == LAB_RAIL) {
            // A single bright leader splits into successively darker return strokes.
            railRumble += .038f * (low - railRumble);
            float fork = t - (.014f + .045f*variation);
            float last = t - (.047f - .08f*variation);
            float snap = 1.65f*expf(-185*t);
            if (fork >= 0) snap += .95f*expf(-125*fork);
            if (last >= 0) snap += .60f*expf(-100*last);
            y = snap*(noise-low) + 1.8f*low*expf(-35*t)
              + 2.3f*railRumble*expf(-17*t);
        } else if (voice == LAB_GLOB) {
            float phase = 2*PI*(105*t + 5.1f*(1-expf(-54*t)));
            y = .59f * (sinf(phase) + .21f*sinf(2*phase)) * expf(-(34+variation*25)*t);
        } else if (voice == LAB_LAUNCH || voice == LAB_BOOM) {
            bool blast = voice == LAB_BOOM;
            railRumble += .035f*(low-railRumble);
            float thud = sinf(2*PI*(54*t + 1.25f*(1-expf(-65*t))));
            y = (.64f*thud + 1.8f*low)*expf(-(blast ? 17 : 50)*t)
              + (blast ? 4.3f*railRumble*expf(-7*t) : .60f*noise*expf(-26*t));
        } else if (voice == LAB_DMR || voice == LAB_GLOCK) {
            bool rifle = voice == LAB_DMR;
            float returnTime = t - (rifle ? .075f : .031f) - variation*.05f;
            y = (rifle ? 1.42f : 1.08f)*(noise-low)*expf(-(rifle ? 104 : 215)*t)
              + (rifle ? 1.9f : .85f)*low*expf(-(rifle ? 42 : 100)*t);
            if (returnTime > 0) {
                float metal = sinf(2*PI*(rifle ? 540 : 1570)*returnTime);
                y += (.38f*noise + .27f*metal)*expf(-(rifle ? 88 : 158)*returnTime);
            }
        } else {
            float phase = 2*PI*(180*t + 2.4f*(1-expf(-57*t)));
            y = (.61f*low + .19f*sinf(phase + .7f*sinf(2*PI*49*t))) * expf(-(27+variation*40)*t);
        }
        float edge = fminf(1,t*2400) * fminf(1,(frames-1-i)*.006f);
        pcm[i] = (short)(30000 * tanhf(y) * edge);
    }
    if (voice == LAB_RAIL) lab_lightning_reverb(pcm, frames, take);
    return (Wave){frames, LAB_RATE, 16, 1, pcm};
}

// Phase and envelope belong to the audio thread. The game publishes only charge.
static void lab_charge_callback(void *buffer, unsigned frames)
{
    float *out = buffer;
    int level = atomic_load(&labChargeLevel);
    float progress = fmaxf(0,(level-1)/1000.f);
    float hz = 244 + 232*progress;
    for (unsigned i=0; i<frames; ++i) {
        labChargeGain += .006f * ((level ? .078f : 0) - labChargeGain);
        labPhaseA += hz/LAB_RATE;
        labPhaseB += (hz + 15*(1-progress))/LAB_RATE;
        labPhaseA -= floorf(labPhaseA);
        labPhaseB -= floorf(labPhaseB);
        // A stronger fundamental keeps the full-charge tone audible at any phase.
        out[i] = labChargeGain * (sinf(2*PI*labPhaseA) + .55f*sinf(2*PI*labPhaseB));
    }
}

// Sustained sources share mono samples. Their positions and motor pitch update in game time.
static void lab_load_sustain(SfxId id)
{
    unsigned frames = LAB_RATE * (id == SFX_LAB_MOTOR ? 6 : 2);
    short *pcm = MemAlloc(frames*sizeof(*pcm));
    if (!pcm) return;
    uint32_t seed = 725;
    float low = 0;
    for (unsigned i=0; i<frames; ++i) {
        float t = (float)i/LAB_RATE, noise = lab_noise(&seed);
        low += .18f*(noise-low);
        float sample;
        if (id == SFX_LAB_BEAM)
            sample = .10f*sinf(2*PI*455*t+.46f*sinf(2*PI*112*t)) + .035f*sinf(2*PI*910*t);
        else if (id == SFX_LAB_SIZZLE)
            sample = (.20f+.11f*sinf(2*PI*39*t))*(noise-low) + .36f*low;
        else
            sample = .61f*low + .12f*(noise-low) + .075f*sinf(2*PI*340*t+.6f*sinf(2*PI*65*t));
        float edge = fminf(1,i/120.f)*fminf(1,(frames-1-i)/120.f);
        pcm[i] = (short)(30000*sample*edge);
    }
    Wave wave = {frames,LAB_RATE,16,1,pcm};
    sfxSounds[id] = LoadSoundFromWave(wave);
    UnloadWave(wave);
}

static void lab_audio_init(void)
{
    for (int voice=0; voice<LAB_VOICE_COUNT; ++voice)
        for (int take=0; take<LAB_VARIANTS; ++take) {
            Wave wave = lab_voice_wave(voice,take);
            sfxSounds[SFX_LAB_FIRST + voice*LAB_VARIANTS + take] = LoadSoundFromWave(wave);
            UnloadWave(wave);
        }
    for (int id=SFX_LAB_BEAM; id<=SFX_LAB_MOTOR; ++id) lab_load_sustain(id);
    labChargeStream = LoadAudioStream(LAB_RATE,32,1);
    SetAudioStreamCallback(labChargeStream,lab_charge_callback);
}
