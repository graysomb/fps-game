#include "net_protocol.h"

#include <assert.h>

int main(void) {
    uint8_t packet[128];
    NetInputCommand input = {
        .sequence = 0x12345678u,
        .client_tick = 91u,
        .move_x = -1234,
        .move_y = 32767,
        .look_x = 222,
        .look_y = -333,
        .held = NET_INPUT_TETHER,
        .pressed = NET_INPUT_FIRE | NET_INPUT_JUMP,
        .creative_brush = 7,
        .creative_color = 4,
        .creative_pickup = 3
    };
    NetWriter writer;
    net_writer_init(&writer, packet, sizeof(packet));
    assert(net_write_header(&writer, NET_MSG_INPUT, 2, 77, 88));
    net_write_u8(&writer, 3);
    assert(net_write_input(&writer, &input));
    assert(!writer.failed);

    NetReader reader;
    NetPacketHeader header;
    NetInputCommand decoded;
    net_reader_init(&reader, packet, writer.length);
    assert(net_read_header(&reader, &header));
    assert(header.type == NET_MSG_INPUT && header.flags == 2);
    assert(header.session_id == 77 && header.server_tick == 88);
    assert(net_read_u8(&reader) == 3);
    assert(net_read_input(&reader, &decoded));
    assert(decoded.sequence == input.sequence && decoded.client_tick == input.client_tick);
    assert(decoded.move_x == input.move_x && decoded.move_y == input.move_y);
    assert(decoded.look_x == input.look_x && decoded.look_y == input.look_y);
    assert(decoded.held == input.held && decoded.pressed == input.pressed);
    assert(decoded.creative_brush == input.creative_brush);
    assert(decoded.creative_color == input.creative_color);
    assert(decoded.creative_pickup == input.creative_pickup);

    NetPlayerVisualState visual = {
        .flags = NET_PLAYER_VISUAL_MELEE | NET_PLAYER_VISUAL_TETHER,
        .melee_progress = 173,
        .tether_x = 1.25f,
        .tether_y = -2.5f,
        .tether_z = 9.75f
    };
    net_writer_init(&writer, packet, sizeof(packet));
    assert(net_write_player_visual(&writer, &visual));
    NetPlayerVisualState decoded_visual;
    net_reader_init(&reader, packet, writer.length);
    assert(net_read_player_visual(&reader, &decoded_visual));
    assert(decoded_visual.flags == visual.flags);
    assert(decoded_visual.melee_progress == visual.melee_progress);
    assert(decoded_visual.tether_x == visual.tether_x);
    assert(decoded_visual.tether_y == visual.tether_y);
    assert(decoded_visual.tether_z == visual.tether_z);

    net_reader_init(&reader, packet, writer.length - 1);
    assert(!net_read_player_visual(&reader, &decoded_visual));

    net_writer_init(&writer, packet, sizeof(packet));
    assert(net_write_header(&writer, NET_MSG_INPUT, 2, 77, 88));
    net_write_u8(&writer, 3);
    assert(net_write_input(&writer, &input));
    net_reader_init(&reader, packet, writer.length - 1);
    assert(net_read_header(&reader, &header));
    assert(net_read_u8(&reader) == 3);
    assert(!net_read_input(&reader, &decoded));

    packet[0] ^= 1u;
    net_reader_init(&reader, packet, writer.length);
    assert(!net_read_header(&reader, &header));
    return 0;
}
