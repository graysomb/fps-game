#ifndef FPS_NET_TRANSPORT_H
#define FPS_NET_TRANSPORT_H

#include "net_protocol.h"

#include <enet/enet.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define FPS_NET_CHANNEL_CONTROL 0u
#define FPS_NET_CHANNEL_WORLD 1u
#define FPS_NET_CHANNEL_INPUT 2u
#define FPS_NET_CHANNEL_SNAPSHOT 3u
#define FPS_NET_CHANNEL_COUNT 4u

typedef enum NetRole {
    NET_ROLE_OFFLINE = 0,
    NET_ROLE_HOST,
    NET_ROLE_CLIENT
} NetRole;

typedef struct NetPeerInfo {
    ENetPeer *peer;
    bool connected;
    bool welcomed;
    bool ready;
    uint32_t last_input_sequence;
    uint32_t last_acked_input_sequence;
    double last_input_time;
} NetPeerInfo;

typedef struct NetTransport {
    NetRole role;
    ENetHost *host;
    ENetPeer *server_peer;
    NetPeerInfo peers[FPS_NET_MAX_PLAYERS];
    uint16_t port;
    uint32_t session_id;
    int local_player_slot;
    bool connected;
    bool initialized;
    char remote_host[256];
    uint64_t bytes_sent;
    uint64_t bytes_received;
    uint64_t packets_sent;
    uint64_t packets_received;
    char last_error[160];
} NetTransport;

typedef void (*NetReceiveFn)(NetTransport *transport, int player_slot,
                             uint8_t channel, const uint8_t *data, size_t length,
                             void *user);
typedef void (*NetConnectFn)(NetTransport *transport, int player_slot,
                             bool connected, void *user);

bool net_transport_global_init(void);
void net_transport_global_shutdown(void);
bool net_transport_host(NetTransport *transport, uint16_t port);
bool net_transport_connect(NetTransport *transport, const char *host, uint16_t port);
void net_transport_pump(NetTransport *transport, NetReceiveFn receive_fn,
                        NetConnectFn connect_fn, void *user);
bool net_transport_send(NetTransport *transport, int player_slot, uint8_t channel,
                        const void *data, size_t length, bool reliable);
void net_transport_broadcast(NetTransport *transport, uint8_t channel,
                             const void *data, size_t length, bool reliable);
void net_transport_disconnect(NetTransport *transport, uint32_t reason);
void net_transport_shutdown(NetTransport *transport);
int net_transport_connected_clients(const NetTransport *transport);
uint32_t net_transport_rtt_ms(const NetTransport *transport);

#endif
