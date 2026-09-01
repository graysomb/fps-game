#include "net_protocol.h"

#include <string.h>

static void net_writer_reserve(NetWriter *writer, size_t amount) {
    if (!writer || writer->failed || amount > writer->capacity - writer->length) {
        if (writer) writer->failed = true;
    }
}

static void net_reader_reserve(NetReader *reader, size_t amount) {
    if (!reader || reader->failed || amount > reader->length - reader->offset) {
        if (reader) reader->failed = true;
    }
}

void net_writer_init(NetWriter *writer, uint8_t *data, size_t capacity) {
    if (!writer) return;
    writer->data = data;
    writer->capacity = capacity;
    writer->length = 0;
    writer->failed = data == NULL;
}

void net_reader_init(NetReader *reader, const uint8_t *data, size_t length) {
    if (!reader) return;
    reader->data = data;
    reader->length = length;
    reader->offset = 0;
    reader->failed = data == NULL;
}

void net_write_u8(NetWriter *writer, uint8_t value) {
    net_writer_reserve(writer, 1);
    if (!writer || writer->failed) return;
    writer->data[writer->length++] = value;
}

void net_write_u16(NetWriter *writer, uint16_t value) {
    net_writer_reserve(writer, 2);
    if (!writer || writer->failed) return;
    writer->data[writer->length++] = (uint8_t)(value >> 8);
    writer->data[writer->length++] = (uint8_t)value;
}

void net_write_u32(NetWriter *writer, uint32_t value) {
    net_writer_reserve(writer, 4);
    if (!writer || writer->failed) return;
    writer->data[writer->length++] = (uint8_t)(value >> 24);
    writer->data[writer->length++] = (uint8_t)(value >> 16);
    writer->data[writer->length++] = (uint8_t)(value >> 8);
    writer->data[writer->length++] = (uint8_t)value;
}

void net_write_u64(NetWriter *writer, uint64_t value) {
    net_write_u32(writer, (uint32_t)(value >> 32));
    net_write_u32(writer, (uint32_t)value);
}

void net_write_i16(NetWriter *writer, int16_t value) { net_write_u16(writer, (uint16_t)value); }
void net_write_i32(NetWriter *writer, int32_t value) { net_write_u32(writer, (uint32_t)value); }

void net_write_f32(NetWriter *writer, float value) {
    uint32_t bits = 0;
    memcpy(&bits, &value, sizeof(bits));
    net_write_u32(writer, bits);
}

void net_write_bytes(NetWriter *writer, const void *data, size_t length) {
    net_writer_reserve(writer, length);
    if (!writer || writer->failed) return;
    memcpy(writer->data + writer->length, data, length);
    writer->length += length;
}

uint8_t net_read_u8(NetReader *reader) {
    net_reader_reserve(reader, 1);
    return (!reader || reader->failed) ? 0 : reader->data[reader->offset++];
}

uint16_t net_read_u16(NetReader *reader) {
    net_reader_reserve(reader, 2);
    if (!reader || reader->failed) return 0;
    uint16_t value = (uint16_t)reader->data[reader->offset] << 8 |
                     (uint16_t)reader->data[reader->offset + 1];
    reader->offset += 2;
    return value;
}

uint32_t net_read_u32(NetReader *reader) {
    net_reader_reserve(reader, 4);
    if (!reader || reader->failed) return 0;
    uint32_t value = (uint32_t)reader->data[reader->offset] << 24 |
                     (uint32_t)reader->data[reader->offset + 1] << 16 |
                     (uint32_t)reader->data[reader->offset + 2] << 8 |
                     (uint32_t)reader->data[reader->offset + 3];
    reader->offset += 4;
    return value;
}

uint64_t net_read_u64(NetReader *reader) {
    uint64_t high = net_read_u32(reader);
    uint64_t low = net_read_u32(reader);
    return (high << 32) | low;
}

int16_t net_read_i16(NetReader *reader) { return (int16_t)net_read_u16(reader); }
int32_t net_read_i32(NetReader *reader) { return (int32_t)net_read_u32(reader); }

float net_read_f32(NetReader *reader) {
    uint32_t bits = net_read_u32(reader);
    float value = 0.0f;
    memcpy(&value, &bits, sizeof(value));
    return value;
}

void net_read_bytes(NetReader *reader, void *data, size_t length) {
    net_reader_reserve(reader, length);
    if (!reader || reader->failed) return;
    memcpy(data, reader->data + reader->offset, length);
    reader->offset += length;
}

bool net_write_header(NetWriter *writer, NetMessageType type, uint8_t flags,
                      uint32_t session_id, uint32_t server_tick) {
    net_write_u32(writer, FPS_NET_MAGIC);
    net_write_u16(writer, FPS_NET_PROTOCOL_VERSION);
    net_write_u8(writer, (uint8_t)type);
    net_write_u8(writer, flags);
    net_write_u32(writer, session_id);
    net_write_u32(writer, server_tick);
    return writer && !writer->failed;
}

bool net_read_header(NetReader *reader, NetPacketHeader *header) {
    if (!reader || !header) return false;
    uint32_t magic = net_read_u32(reader);
    uint16_t version = net_read_u16(reader);
    header->type = net_read_u8(reader);
    header->flags = net_read_u8(reader);
    header->session_id = net_read_u32(reader);
    header->server_tick = net_read_u32(reader);
    return !reader->failed && magic == FPS_NET_MAGIC && version == FPS_NET_PROTOCOL_VERSION;
}

bool net_write_input(NetWriter *writer, const NetInputCommand *command) {
    if (!writer || !command) return false;
    net_write_u32(writer, command->sequence);
    net_write_u32(writer, command->client_tick);
    net_write_i16(writer, command->move_x);
    net_write_i16(writer, command->move_y);
    net_write_i16(writer, command->look_x);
    net_write_i16(writer, command->look_y);
    net_write_u16(writer, command->held);
    net_write_u16(writer, command->pressed);
    net_write_u8(writer, command->creative_brush);
    net_write_u8(writer, command->creative_color);
    net_write_u8(writer, command->creative_pickup);
    return !writer->failed;
}

bool net_read_input(NetReader *reader, NetInputCommand *command) {
    if (!reader || !command) return false;
    memset(command, 0, sizeof(*command));
    command->sequence = net_read_u32(reader);
    command->client_tick = net_read_u32(reader);
    command->move_x = net_read_i16(reader);
    command->move_y = net_read_i16(reader);
    command->look_x = net_read_i16(reader);
    command->look_y = net_read_i16(reader);
    command->held = net_read_u16(reader);
    command->pressed = net_read_u16(reader);
    command->creative_brush = net_read_u8(reader);
    command->creative_color = net_read_u8(reader);
    command->creative_pickup = net_read_u8(reader);
    return !reader->failed;
}
