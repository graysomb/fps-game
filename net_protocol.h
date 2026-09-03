#ifndef FPS_NET_PROTOCOL_H
#define FPS_NET_PROTOCOL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define FPS_NET_MAGIC 0x4650534eu /* "FPSN" */
#define FPS_NET_PROTOCOL_VERSION 3u
#define FPS_NET_DEFAULT_PORT 27015u
#define FPS_NET_DISCOVERY_PORT 27016u
#define FPS_NET_MAX_PACKET 1100u
#define FPS_NET_TICK_RATE 60u
#define FPS_NET_POSE_RATE 60u
#define FPS_NET_MAX_PLAYERS 4u

typedef enum NetMessageType {
    NET_MSG_DISCOVER = 1,
    NET_MSG_DISCOVER_REPLY,
    NET_MSG_HELLO,
    NET_MSG_WELCOME,
    NET_MSG_REJECT,
    NET_MSG_INPUT,
    NET_MSG_PLAYER_STATE,
    NET_MSG_WORLD_BEGIN,
    NET_MSG_WORLD_STATIC,
    NET_MSG_WORLD_PICKUPS,
    NET_MSG_WORLD_END,
    NET_MSG_READY,
    NET_MSG_START,
    NET_MSG_STATIC_RESET,
    NET_MSG_DYNAMIC_POSES,
    NET_MSG_DYNAMIC_DELETE,
    NET_MSG_SESSION_STATE,
    NET_MSG_CREATIVE_ACTION,
    NET_MSG_GAME_EVENT
} NetMessageType;

typedef enum NetInputButtons {
    NET_INPUT_JUMP = 1u << 0,
    NET_INPUT_FIRE = 1u << 1,
    NET_INPUT_MELEE = 1u << 2,
    NET_INPUT_BUILD = 1u << 3,
    NET_INPUT_TETHER = 1u << 4,
    NET_INPUT_CREATIVE_PLACE = 1u << 5,
    NET_INPUT_CREATIVE_REMOVE = 1u << 6,
    NET_INPUT_CREATIVE_PICKUP = 1u << 7,
    NET_INPUT_CREATIVE_REMOVE_PICKUP = 1u << 8,
    NET_INPUT_CREATIVE_ACTIVATE = 1u << 9,
    NET_INPUT_CREATIVE_UP = 1u << 10,
    NET_INPUT_CREATIVE_DOWN = 1u << 11
} NetInputButtons;

typedef struct NetPacketHeader {
    uint8_t type;
    uint8_t flags;
    uint32_t session_id;
    uint32_t server_tick;
} NetPacketHeader;

typedef struct NetInputCommand {
    uint32_t sequence;
    uint32_t client_tick;
    int16_t move_x;
    int16_t move_y;
    int16_t look_x;
    int16_t look_y;
    uint16_t held;
    uint16_t pressed;
    uint8_t creative_brush;
    uint8_t creative_color;
    uint8_t creative_pickup;
} NetInputCommand;

typedef enum NetPlayerVisualFlags {
    NET_PLAYER_VISUAL_MELEE = 1u << 0,
    NET_PLAYER_VISUAL_TETHER = 1u << 1
} NetPlayerVisualFlags;

typedef struct NetPlayerVisualState {
    uint8_t flags;
    uint8_t melee_progress;
    float tether_x;
    float tether_y;
    float tether_z;
} NetPlayerVisualState;

typedef struct NetWriter {
    uint8_t *data;
    size_t capacity;
    size_t length;
    bool failed;
} NetWriter;

typedef struct NetReader {
    const uint8_t *data;
    size_t length;
    size_t offset;
    bool failed;
} NetReader;

void net_writer_init(NetWriter *writer, uint8_t *data, size_t capacity);
void net_reader_init(NetReader *reader, const uint8_t *data, size_t length);
bool net_write_header(NetWriter *writer, NetMessageType type, uint8_t flags,
                      uint32_t session_id, uint32_t server_tick);
bool net_read_header(NetReader *reader, NetPacketHeader *header);

void net_write_u8(NetWriter *writer, uint8_t value);
void net_write_u16(NetWriter *writer, uint16_t value);
void net_write_u32(NetWriter *writer, uint32_t value);
void net_write_u64(NetWriter *writer, uint64_t value);
void net_write_i16(NetWriter *writer, int16_t value);
void net_write_i32(NetWriter *writer, int32_t value);
void net_write_f32(NetWriter *writer, float value);
void net_write_bytes(NetWriter *writer, const void *data, size_t length);

uint8_t net_read_u8(NetReader *reader);
uint16_t net_read_u16(NetReader *reader);
uint32_t net_read_u32(NetReader *reader);
uint64_t net_read_u64(NetReader *reader);
int16_t net_read_i16(NetReader *reader);
int32_t net_read_i32(NetReader *reader);
float net_read_f32(NetReader *reader);
void net_read_bytes(NetReader *reader, void *data, size_t length);

bool net_write_input(NetWriter *writer, const NetInputCommand *command);
bool net_read_input(NetReader *reader, NetInputCommand *command);
bool net_write_player_visual(NetWriter *writer, const NetPlayerVisualState *state);
bool net_read_player_visual(NetReader *reader, NetPlayerVisualState *state);

#endif
