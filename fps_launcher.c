#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FPS_EXIT_GPU_CONTEXT_UNAVAILABLE 78
#define FPS_EXIT_GPU_INITIALIZATION_FAILED 79

#ifdef _WIN32
#include <windows.h>
#include <process.h>
#include <direct.h>
#define PATH_SEPARATOR '\\'
#define GPU_BINARY "fps_ray_gpu.exe"
#define CPU_BINARY "fps_ray_cpu.exe"
#else
#include <limits.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#define PATH_SEPARATOR '/'
#define GPU_BINARY "fps_ray_gpu"
#define CPU_BINARY "fps_ray_cpu"
#endif

static bool requests_cpu(int argc, char **argv) {
    for (int i = 1; i < argc; ++i) {
        const char *value = NULL;
        if (strncmp(argv[i], "--physics=", 10) == 0) value = argv[i] + 10;
        else if (strcmp(argv[i], "--physics") == 0 && i + 1 < argc) value = argv[++i];
        if (value && (strcmp(value, "cpu-mt") == 0 || strcmp(value, "cpu-st") == 0 ||
                      strcmp(value, "mt") == 0 || strcmp(value, "st") == 0)) return true;
    }
    return false;
}

static bool forces_gpu(int argc, char **argv) {
    for (int i = 1; i < argc; ++i) {
        const char *value = NULL;
        if (strncmp(argv[i], "--physics=", 10) == 0) value = argv[i] + 10;
        else if (strcmp(argv[i], "--physics") == 0 && i + 1 < argc) value = argv[++i];
        if (value && (strcmp(value, "gpu") == 0 || strcmp(value, "gpu-gl43") == 0)) return true;
    }
    return false;
}

static bool executable_directory(char *out, size_t capacity) {
#ifdef _WIN32
    DWORD length = GetModuleFileNameA(NULL, out, (DWORD)capacity);
    if (length == 0 || length >= capacity) return false;
#else
    ssize_t length = readlink("/proc/self/exe", out, capacity - 1);
    if (length <= 0 || (size_t)length >= capacity) return false;
    out[length] = '\0';
#endif
    char *slash = strrchr(out, PATH_SEPARATOR);
    if (!slash) return false;
    *slash = '\0';
    return true;
}

static int run_child(const char *directory, const char *binary, int argc, char **argv) {
    size_t path_len = strlen(directory) + strlen(binary) + 2;
    char *path = malloc(path_len);
    char **child_argv = calloc((size_t)argc + 1, sizeof(*child_argv));
    if (!path || !child_argv) {
        free(path); free(child_argv);
        return -1;
    }
    snprintf(path, path_len, "%s%c%s", directory, PATH_SEPARATOR, binary);
    child_argv[0] = path;
    for (int i = 1; i < argc; ++i) child_argv[i] = argv[i];
    child_argv[argc] = NULL;
#ifdef _WIN32
    intptr_t result = _spawnv(_P_WAIT, path, (const char *const *)child_argv);
    int status = (result < 0) ? -1 : (int)result;
#else
    pid_t pid = fork();
    if (pid == 0) {
        execv(path, child_argv);
        _exit(127);
    }
    int status = -1;
    if (pid > 0) {
        int wait_status = 0;
        if (waitpid(pid, &wait_status, 0) >= 0) {
            if (WIFEXITED(wait_status)) status = WEXITSTATUS(wait_status);
            else if (WIFSIGNALED(wait_status)) status = 128 + WTERMSIG(wait_status);
        }
    }
#endif
    free(path); free(child_argv);
    return status;
}

int main(int argc, char **argv) {
    char directory[4096];
    if (!executable_directory(directory, sizeof(directory))) {
        fprintf(stderr, "Unable to locate FPS game binaries\n");
        return 1;
    }
#ifdef _WIN32
    if (_chdir(directory) != 0) {
#else
    if (chdir(directory) != 0) {
#endif
        fprintf(stderr, "Unable to use executable directory as working directory: %s\n", strerror(errno));
        return 1;
    }
    if (requests_cpu(argc, argv)) {
        int status = run_child(directory, CPU_BINARY, argc, argv);
        if (status < 0) fprintf(stderr, "Unable to start %s: %s\n", CPU_BINARY, strerror(errno));
        return (status < 0) ? 1 : status;
    }
    int status = run_child(directory, GPU_BINARY, argc, argv);
    if (forces_gpu(argc, argv)) {
        if (status < 0) fprintf(stderr, "Unable to start %s: %s\n", GPU_BINARY, strerror(errno));
        return (status < 0) ? 1 : status;
    }
    if (status < 0 || status == FPS_EXIT_GPU_CONTEXT_UNAVAILABLE || status == FPS_EXIT_GPU_INITIALIZATION_FAILED) {
        fprintf(stderr, "GPU build unavailable; starting CPU compatibility build\n");
        status = run_child(directory, CPU_BINARY, argc, argv);
    }
    if (status < 0) fprintf(stderr, "Unable to start %s: %s\n", CPU_BINARY, strerror(errno));
    return (status < 0) ? 1 : status;
}
