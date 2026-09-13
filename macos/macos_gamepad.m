#define MACOS_GAMEPAD_IMPLEMENTATION
#include "macos/macos_gamepad.h"

#include "raylib.h"

#import <Foundation/Foundation.h>
#import <GameController/GameController.h>
#import <CoreFoundation/CoreFoundation.h>

#include <string.h>

int glfwJoystickIsGamepad(int jid);

#define MACOS_GAMEPAD_MAX 4
#define MACOS_GAMEPAD_BUTTONS 18
#define MACOS_GAMEPAD_AXES 6

static bool gamepad_inited;
static bool have_snapshot;
static bool using_gc;
static bool present[MACOS_GAMEPAD_MAX];
static bool button_down[MACOS_GAMEPAD_MAX][MACOS_GAMEPAD_BUTTONS];
static bool button_prev[MACOS_GAMEPAD_MAX][MACOS_GAMEPAD_BUTTONS];
static float axes[MACOS_GAMEPAD_MAX][MACOS_GAMEPAD_AXES];

static void reset_axes(void) {
    memset(axes, 0, sizeof(axes));
    for (int i = 0; i < MACOS_GAMEPAD_MAX; i++) {
        axes[i][GAMEPAD_AXIS_LEFT_TRIGGER] = -1.0f;
        axes[i][GAMEPAD_AXIS_RIGHT_TRIGGER] = -1.0f;
    }
}

static float analog_to_glfw_trigger(float value) {
    if (value < 0.0f) value = 0.0f;
    if (value > 1.0f) value = 1.0f;
    return value * 2.0f - 1.0f;
}

static void snapshot_from_gc(void) {
    memcpy(button_prev, button_down, sizeof(button_down));
    memset(button_down, 0, sizeof(button_down));
    memset(present, 0, sizeof(present));
    reset_axes();
    using_gc = false;

    NSArray<GCController *> *controllers = [GCController controllers];
    if (controllers.count == 0) return;

    int slot = 0;
    for (GCController *controller in controllers) {
        if (slot >= MACOS_GAMEPAD_MAX) break;
        GCExtendedGamepad *pad = controller.extendedGamepad;
        if (!pad) continue;

        using_gc = true;
        present[slot] = true;
        if (controller.playerIndex == GCControllerPlayerIndexUnset) {
            controller.playerIndex = (GCControllerPlayerIndex)slot;
        }
        if (!pad.valueChangedHandler) {
            pad.valueChangedHandler = ^(GCExtendedGamepad *gp, GCControllerElement *element) {
                (void)gp;
                (void)element;
            };
        }

        // GameController stick Y is up-positive; GLFW/raylib is down-positive.
        axes[slot][GAMEPAD_AXIS_LEFT_X] = pad.leftThumbstick.xAxis.value;
        axes[slot][GAMEPAD_AXIS_LEFT_Y] = -pad.leftThumbstick.yAxis.value;
        axes[slot][GAMEPAD_AXIS_RIGHT_X] = pad.rightThumbstick.xAxis.value;
        axes[slot][GAMEPAD_AXIS_RIGHT_Y] = -pad.rightThumbstick.yAxis.value;
        axes[slot][GAMEPAD_AXIS_LEFT_TRIGGER] = analog_to_glfw_trigger(pad.leftTrigger.value);
        axes[slot][GAMEPAD_AXIS_RIGHT_TRIGGER] = analog_to_glfw_trigger(pad.rightTrigger.value);

        button_down[slot][GAMEPAD_BUTTON_RIGHT_FACE_DOWN] = pad.buttonA.pressed;
        button_down[slot][GAMEPAD_BUTTON_RIGHT_FACE_RIGHT] = pad.buttonB.pressed;
        button_down[slot][GAMEPAD_BUTTON_RIGHT_FACE_LEFT] = pad.buttonX.pressed;
        button_down[slot][GAMEPAD_BUTTON_RIGHT_FACE_UP] = pad.buttonY.pressed;
        button_down[slot][GAMEPAD_BUTTON_LEFT_TRIGGER_1] = pad.leftShoulder.pressed;
        button_down[slot][GAMEPAD_BUTTON_RIGHT_TRIGGER_1] = pad.rightShoulder.pressed;
        button_down[slot][GAMEPAD_BUTTON_LEFT_TRIGGER_2] = pad.leftTrigger.pressed;
        button_down[slot][GAMEPAD_BUTTON_RIGHT_TRIGGER_2] = pad.rightTrigger.pressed;
        button_down[slot][GAMEPAD_BUTTON_LEFT_FACE_UP] = pad.dpad.up.pressed;
        button_down[slot][GAMEPAD_BUTTON_LEFT_FACE_RIGHT] = pad.dpad.right.pressed;
        button_down[slot][GAMEPAD_BUTTON_LEFT_FACE_DOWN] = pad.dpad.down.pressed;
        button_down[slot][GAMEPAD_BUTTON_LEFT_FACE_LEFT] = pad.dpad.left.pressed;
        button_down[slot][GAMEPAD_BUTTON_LEFT_THUMB] = pad.leftThumbstickButton.pressed;
        button_down[slot][GAMEPAD_BUTTON_RIGHT_THUMB] = pad.rightThumbstickButton.pressed;
        button_down[slot][GAMEPAD_BUTTON_MIDDLE_RIGHT] = pad.buttonMenu.pressed;
        button_down[slot][GAMEPAD_BUTTON_MIDDLE_LEFT] = pad.buttonOptions.pressed;
        button_down[slot][GAMEPAD_BUTTON_MIDDLE] = pad.buttonHome.pressed;
        slot++;
    }
}

void macos_gamepad_init(void) {
    if (gamepad_inited) return;
    gamepad_inited = true;
    reset_axes();
    [GCController startWirelessControllerDiscoveryWithCompletionHandler:nil];
    CFRunLoopRunInMode(kCFRunLoopDefaultMode, 0.0, false);
}

void macos_gamepad_poll(void) {
    macos_gamepad_init();
    @autoreleasepool {
        snapshot_from_gc();
    }
    have_snapshot = true;
}

static void ensure_poll(void) {
    if (!have_snapshot) macos_gamepad_poll();
}

static bool valid_pad(int gamepad) {
    return gamepad >= 0 && gamepad < MACOS_GAMEPAD_MAX && present[gamepad];
}

static bool raylib_mapped_pad(int gamepad) {
    return IsGamepadAvailable(gamepad) && glfwJoystickIsGamepad(gamepad);
}

bool macos_IsGamepadAvailable(int gamepad) {
    ensure_poll();
    if (using_gc) return valid_pad(gamepad);
    return raylib_mapped_pad(gamepad);
}

float macos_GetGamepadAxisMovement(int gamepad, int axis) {
    ensure_poll();
    if (using_gc) {
        if (!valid_pad(gamepad) || axis < 0 || axis >= MACOS_GAMEPAD_AXES) {
            return (axis == GAMEPAD_AXIS_LEFT_TRIGGER || axis == GAMEPAD_AXIS_RIGHT_TRIGGER) ? -1.0f : 0.0f;
        }
        return axes[gamepad][axis];
    }
    if (!raylib_mapped_pad(gamepad)) {
        return (axis == GAMEPAD_AXIS_LEFT_TRIGGER || axis == GAMEPAD_AXIS_RIGHT_TRIGGER) ? -1.0f : 0.0f;
    }
    return GetGamepadAxisMovement(gamepad, axis);
}

bool macos_IsGamepadButtonDown(int gamepad, int button) {
    ensure_poll();
    if (using_gc) {
        if (!valid_pad(gamepad) || button <= 0 || button >= MACOS_GAMEPAD_BUTTONS) return false;
        return button_down[gamepad][button];
    }
    return raylib_mapped_pad(gamepad) && IsGamepadButtonDown(gamepad, button);
}

bool macos_IsGamepadButtonUp(int gamepad, int button) {
    ensure_poll();
    if (using_gc) {
        if (!valid_pad(gamepad) || button <= 0 || button >= MACOS_GAMEPAD_BUTTONS) return true;
        return !button_down[gamepad][button];
    }
    if (!raylib_mapped_pad(gamepad)) return true;
    return IsGamepadButtonUp(gamepad, button);
}

bool macos_IsGamepadButtonPressed(int gamepad, int button) {
    ensure_poll();
    if (using_gc) {
        if (!valid_pad(gamepad) || button <= 0 || button >= MACOS_GAMEPAD_BUTTONS) return false;
        return button_down[gamepad][button] && !button_prev[gamepad][button];
    }
    return raylib_mapped_pad(gamepad) && IsGamepadButtonPressed(gamepad, button);
}

bool macos_IsGamepadButtonReleased(int gamepad, int button) {
    ensure_poll();
    if (using_gc) {
        if (!valid_pad(gamepad) || button <= 0 || button >= MACOS_GAMEPAD_BUTTONS) return false;
        return !button_down[gamepad][button] && button_prev[gamepad][button];
    }
    return raylib_mapped_pad(gamepad) && IsGamepadButtonReleased(gamepad, button);
}
