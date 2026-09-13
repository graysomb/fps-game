#ifndef MACOS_GAMEPAD_H
#define MACOS_GAMEPAD_H

#include <stdbool.h>

void macos_gamepad_init(void);
void macos_gamepad_poll(void);
bool macos_IsGamepadAvailable(int gamepad);
float macos_GetGamepadAxisMovement(int gamepad, int axis);
bool macos_IsGamepadButtonPressed(int gamepad, int button);
bool macos_IsGamepadButtonDown(int gamepad, int button);
bool macos_IsGamepadButtonReleased(int gamepad, int button);
bool macos_IsGamepadButtonUp(int gamepad, int button);

#if !defined(MACOS_GAMEPAD_IMPLEMENTATION)
#define IsGamepadAvailable(gamepad) macos_IsGamepadAvailable(gamepad)
#define GetGamepadAxisMovement(gamepad, axis) macos_GetGamepadAxisMovement((gamepad), (axis))
#define IsGamepadButtonPressed(gamepad, button) macos_IsGamepadButtonPressed((gamepad), (button))
#define IsGamepadButtonDown(gamepad, button) macos_IsGamepadButtonDown((gamepad), (button))
#define IsGamepadButtonReleased(gamepad, button) macos_IsGamepadButtonReleased((gamepad), (button))
#define IsGamepadButtonUp(gamepad, button) macos_IsGamepadButtonUp((gamepad), (button))
#endif

#endif
