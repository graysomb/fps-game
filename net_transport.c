#include "net_transport.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int net_global_refs = 0;

static void set_error(NetTransport *transport, const char *message) {
    if (!transport) return;
    snprintf(transport->last_error, sizeof(transport->last_error), "%s",
             message ? message : "network error");
}

bool net_transport_global_init(void) {
    if (net_global_refs++ > 0) return true;
    if (enet_initialize() != 0) {
        net_global_refs = 0;
        return false;
    }
    return true;
}

void net_transport_global_shutdown(void) {
    if (net_global_refs <= 0) return;
    if (--net_global_refs == 0) enet_deinitialize();
}

static void net_transport_reset(NetTransport *transport) {
    if (!transport) return;
    memset(transport, 0, sizeof(*transport));
    transport->role = NET_ROLE_OFFLINE;
    transport->local_player_slot = -1;
}

bool net_transport_host(NetTransport *transport, uint16_t port) {
    if (!transport || !net_transport_global_init()) return false;
    net_transport_reset(transport);
    ENetAddress address;
    memset(&address, 0, sizeof(address));
    address.host = ENET_HOST_ANY;
    address.port = port;
    transport->host = enet_host_create(&address, FPS_NET_MAX_PLAYERS - 1,
                                       FPS_NET_CHANNEL_COUNT, 0, 0);
    if (!transport->host) {
        set_error(transport, "unable to bind LAN host port");
        net_transport_global_shutdown();
        return false;
    }
    transport->role = NET_ROLE_HOST;
    transport->port = port;
    transport->local_player_slot = 0;
    transport->connected = true;
    transport->initialized = true;
    transport->session_id = (uint32_t)enet_time_get() ^ 0x9e3779b9u;
    if (transport->session_id == 0) transport->session_id = 1;
    return true;
}

bool net_transport_connect(NetTransport *transport, const char *host_name, uint16_t port) {
    if (!transport || !host_name || !host_name[0] || !net_transport_global_init()) return false;
    net_transport_reset(transport);
    transport->host = enet_host_create(NULL, 1, FPS_NET_CHANNEL_COUNT, 0, 0);
    if (!transport->host) {
        set_error(transport, "unable to create LAN client");
        net_transport_global_shutdown();
        return false;
    }
    ENetAddress address;
    memset(&address, 0, sizeof(address));
    address.port = port;
    if (enet_address_set_host(&address, host_name) != 0) {
        set_error(transport, "unable to resolve host address");
        enet_host_destroy(transport->host);
        transport->host = NULL;
        net_transport_global_shutdown();
        return false;
    }
    transport->server_peer = enet_host_connect(transport->host, &address,
                                               FPS_NET_CHANNEL_COUNT, 0);
    if (!transport->server_peer) {
        set_error(transport, "unable to start connection");
        enet_host_destroy(transport->host);
        transport->host = NULL;
        net_transport_global_shutdown();
        return false;
    }
    transport->role = NET_ROLE_CLIENT;
    transport->port = port;
    transport->local_player_slot = -1;
    transport->initialized = true;
    snprintf(transport->remote_host, sizeof(transport->remote_host), "%s", host_name);
    return true;
}

static int allocate_host_slot(NetTransport *transport, ENetPeer *peer) {
    for (int slot = 1; slot < (int)FPS_NET_MAX_PLAYERS; ++slot) {
        if (transport->peers[slot].connected) continue;
        transport->peers[slot].peer = peer;
        transport->peers[slot].connected = true;
        transport->peers[slot].last_input_time = (double)enet_time_get() / 1000.0;
        peer->data = (void *)(intptr_t)(slot + 1);
        return slot;
    }
    return -1;
}

static int peer_slot(const NetTransport *transport, const ENetPeer *peer) {
    if (!transport || !peer) return -1;
    if (transport->role == NET_ROLE_CLIENT) return 0;
    intptr_t encoded = (intptr_t)peer->data;
    int slot = (int)encoded - 1;
    return slot >= 1 && slot < (int)FPS_NET_MAX_PLAYERS ? slot : -1;
}

void net_transport_pump(NetTransport *transport, NetReceiveFn receive_fn,
                        NetConnectFn connect_fn, void *user) {
    if (!transport || !transport->host) return;
    ENetEvent event;
    while (enet_host_service(transport->host, &event, 0) > 0) {
        int slot = peer_slot(transport, event.peer);
        switch (event.type) {
            case ENET_EVENT_TYPE_CONNECT:
                if (transport->role == NET_ROLE_HOST) {
                    slot = allocate_host_slot(transport, event.peer);
                    if (slot < 0) {
                        enet_peer_disconnect_now(event.peer, 1);
                        break;
                    }
                } else {
                    transport->connected = true;
                    transport->server_peer = event.peer;
                    slot = 0;
                }
                if (connect_fn) connect_fn(transport, slot, true, user);
                break;
            case ENET_EVENT_TYPE_RECEIVE:
                transport->bytes_received += event.packet->dataLength;
                transport->packets_received++;
                if (receive_fn && event.packet->dataLength <= FPS_NET_MAX_PACKET) {
                    receive_fn(transport, slot, event.channelID,
                               event.packet->data, event.packet->dataLength, user);
                }
                enet_packet_destroy(event.packet);
                break;
            case ENET_EVENT_TYPE_DISCONNECT:
                if (transport->role == NET_ROLE_HOST && slot >= 1) {
                    memset(&transport->peers[slot], 0, sizeof(transport->peers[slot]));
                } else if (transport->role == NET_ROLE_CLIENT) {
                    transport->connected = false;
                    transport->server_peer = NULL;
                    transport->local_player_slot = -1;
                }
                event.peer->data = NULL;
                if (connect_fn) connect_fn(transport, slot, false, user);
                break;
            default:
                break;
        }
    }
}

bool net_transport_send(NetTransport *transport, int player_slot, uint8_t channel,
                        const void *data, size_t length, bool reliable) {
    if (!transport || !transport->host || !data || length == 0 ||
        length > FPS_NET_MAX_PACKET || channel >= FPS_NET_CHANNEL_COUNT) return false;
    ENetPeer *peer = NULL;
    if (transport->role == NET_ROLE_CLIENT) peer = transport->server_peer;
    else if (player_slot >= 1 && player_slot < (int)FPS_NET_MAX_PLAYERS)
        peer = transport->peers[player_slot].peer;
    if (!peer || peer->state != ENET_PEER_STATE_CONNECTED) return false;
    enet_uint32 flags = reliable ? ENET_PACKET_FLAG_RELIABLE : 0;
    ENetPacket *packet = enet_packet_create(data, length, flags);
    if (!packet || enet_peer_send(peer, channel, packet) != 0) {
        if (packet) enet_packet_destroy(packet);
        return false;
    }
    transport->bytes_sent += length;
    transport->packets_sent++;
    return true;
}

void net_transport_broadcast(NetTransport *transport, uint8_t channel,
                             const void *data, size_t length, bool reliable) {
    if (!transport || transport->role != NET_ROLE_HOST) return;
    for (int slot = 1; slot < (int)FPS_NET_MAX_PLAYERS; ++slot) {
        if (transport->peers[slot].connected)
            net_transport_send(transport, slot, channel, data, length, reliable);
    }
}

void net_transport_disconnect(NetTransport *transport, uint32_t reason) {
    if (!transport || !transport->host) return;
    if (transport->role == NET_ROLE_CLIENT && transport->server_peer) {
        enet_peer_disconnect(transport->server_peer, reason);
    } else if (transport->role == NET_ROLE_HOST) {
        for (int slot = 1; slot < (int)FPS_NET_MAX_PLAYERS; ++slot) {
            if (transport->peers[slot].peer) enet_peer_disconnect(transport->peers[slot].peer, reason);
        }
    }
    enet_host_flush(transport->host);
}

void net_transport_shutdown(NetTransport *transport) {
    if (!transport || !transport->initialized) return;
    net_transport_disconnect(transport, 0);
    if (transport->host) enet_host_destroy(transport->host);
    net_transport_reset(transport);
    net_transport_global_shutdown();
}

int net_transport_connected_clients(const NetTransport *transport) {
    if (!transport || transport->role != NET_ROLE_HOST) return 0;
    int count = 0;
    for (int slot = 1; slot < (int)FPS_NET_MAX_PLAYERS; ++slot)
        if (transport->peers[slot].connected) ++count;
    return count;
}

uint32_t net_transport_rtt_ms(const NetTransport *transport) {
    if (!transport) return 0;
    if (transport->role == NET_ROLE_CLIENT && transport->server_peer)
        return transport->server_peer->roundTripTime;
    uint32_t total = 0, count = 0;
    for (int slot = 1; slot < (int)FPS_NET_MAX_PLAYERS; ++slot) {
        if (!transport->peers[slot].connected || !transport->peers[slot].peer) continue;
        total += transport->peers[slot].peer->roundTripTime;
        ++count;
    }
    return count ? total / count : 0;
}
