#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
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
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#define PATH_SEPARATOR '/'
#define GPU_BINARY "fps_ray_gpu"
#define CPU_BINARY "fps_ray_cpu"
#endif

#define MATRIX_PATH_CAPACITY 4096

typedef struct MatrixResult {
    const char *backend;
    const char *binary;
    int exit_status;
    char output_dir[MATRIX_PATH_CAPACITY];
    char report_path[MATRIX_PATH_CAPACITY];
} MatrixResult;

static bool join_path(char *out, size_t capacity, const char *left, const char *right) {
    size_t left_len = strlen(left);
    bool separator = left_len > 0 && left[left_len - 1] != '/' && left[left_len - 1] != '\\';
    int written = snprintf(out, capacity, "%s%s%s", left, separator ? "/" : "", right);
    return written >= 0 && (size_t)written < capacity;
}

static bool make_directory_recursive(const char *path) {
    if (!path || !path[0] || strlen(path) >= MATRIX_PATH_CAPACITY) return false;
    char copy[MATRIX_PATH_CAPACITY];
    strcpy(copy, path);
    for (char *cursor = copy + 1; *cursor; ++cursor) {
        if (*cursor != '/' && *cursor != '\\') continue;
#ifdef _WIN32
        if (cursor == copy + 2 && copy[1] == ':') continue;
#endif
        char saved = *cursor;
        *cursor = '\0';
#ifdef _WIN32
        if (_mkdir(copy) != 0 && errno != EEXIST) return false;
#else
        if (mkdir(copy, 0755) != 0 && errno != EEXIST) return false;
#endif
        *cursor = saved;
    }
#ifdef _WIN32
    return _mkdir(copy) == 0 || errno == EEXIST;
#else
    return mkdir(copy, 0755) == 0 || errno == EEXIST;
#endif
}

static void write_json_string(FILE *file, const char *value) {
    fputc('"', file);
    const unsigned char *cursor = (const unsigned char *)(value ? value : "");
    while (*cursor) {
        unsigned char ch = *cursor++;
        if (ch == '"') fputs("\\\"", file);
        else if (ch == '\\') fputs("\\\\", file);
        else if (ch == '\n') fputs("\\n", file);
        else if (ch == '\r') fputs("\\r", file);
        else if (ch == '\t') fputs("\\t", file);
        else if (ch < 0x20) fprintf(file, "\\u%04x", (unsigned int)ch);
        else fputc((int)ch, file);
    }
    fputc('"', file);
}

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

#ifdef _WIN32
static char *quote_windows_argument(const char *argument) {
    if (!argument) return NULL;
    bool needs_quotes = argument[0] == '\0' || strpbrk(argument, " \t\n\v\"") != NULL;
    if (!needs_quotes) {
        size_t length = strlen(argument) + 1;
        char *copy = malloc(length);
        if (copy) memcpy(copy, argument, length);
        return copy;
    }
    size_t length = strlen(argument);
    if (length > (SIZE_MAX - 3) / 2) return NULL;
    char *quoted = malloc(length * 2 + 3);
    if (!quoted) return NULL;
    char *out = quoted;
    *out++ = '"';
    const char *cursor = argument;
    while (*cursor) {
        size_t slashes = 0;
        while (*cursor == '\\') {
            ++slashes;
            ++cursor;
        }
        if (*cursor == '"') {
            for (size_t i = 0; i < slashes * 2 + 1; ++i) *out++ = '\\';
            *out++ = *cursor++;
        } else {
            size_t copies = *cursor ? slashes : slashes * 2;
            for (size_t i = 0; i < copies; ++i) *out++ = '\\';
            if (*cursor) *out++ = *cursor++;
        }
    }
    *out++ = '"';
    *out = '\0';
    return quoted;
}
#endif

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
#ifdef _WIN32
    char **owned_argv = calloc((size_t)argc, sizeof(*owned_argv));
    if (!owned_argv) {
        free(path); free(child_argv);
        return -1;
    }
    owned_argv[0] = quote_windows_argument(path);
    child_argv[0] = owned_argv[0];
    for (int i = 1; i < argc; ++i) {
        owned_argv[i] = quote_windows_argument(argv[i]);
        child_argv[i] = owned_argv[i];
    }
    child_argv[argc] = NULL;
    for (int i = 0; i < argc; ++i) {
        if (!owned_argv[i]) {
            for (int j = 0; j < argc; ++j) free(owned_argv[j]);
            free(owned_argv); free(path); free(child_argv);
            return -1;
        }
    }
    intptr_t result = _spawnv(_P_WAIT, path, (const char *const *)child_argv);
    int status = (result < 0) ? -1 : (int)result;
    for (int i = 0; i < argc; ++i) free(owned_argv[i]);
    free(owned_argv);
#else
    child_argv[0] = path;
    for (int i = 1; i < argc; ++i) child_argv[i] = argv[i];
    child_argv[argc] = NULL;
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

static bool matrix_requested(int argc, char **argv) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--debug-matrix") == 0) return true;
    }
    return false;
}

static const char *option_value(int argc, char **argv, const char *equals_prefix,
                                const char *separate) {
    size_t prefix_len = strlen(equals_prefix);
    for (int i = 1; i < argc; ++i) {
        if (strncmp(argv[i], equals_prefix, prefix_len) == 0) return argv[i] + prefix_len;
        if (strcmp(argv[i], separate) == 0 && i + 1 < argc) return argv[i + 1];
    }
    return NULL;
}

static bool matrix_argument_is_filtered(int argc, char **argv, int *index) {
    const char *argument = argv[*index];
    if (strcmp(argument, "--debug-matrix") == 0) return true;
    if (strncmp(argument, "--physics=", 10) == 0 ||
        strncmp(argument, "--debug-output=", 15) == 0) return true;
    if (strcmp(argument, "--physics") == 0 || strcmp(argument, "--debug-output") == 0) {
        if (*index + 1 < argc) ++(*index);
        return true;
    }
    return false;
}

static int run_matrix_child(const char *directory, int argc, char **argv,
                            MatrixResult *result) {
    char physics_argument[64];
    char output_argument[MATRIX_PATH_CAPACITY + 32];
    snprintf(physics_argument, sizeof(physics_argument), "--physics=%s", result->backend);
    int output_written = snprintf(output_argument, sizeof(output_argument),
                                  "--debug-output=%s", result->output_dir);
    if (output_written < 0 || (size_t)output_written >= sizeof(output_argument)) return -1;

    char **child_args = calloc((size_t)argc + 3, sizeof(*child_args));
    if (!child_args) return -1;
    int child_argc = 0;
    child_args[child_argc++] = argv[0];
    for (int i = 1; i < argc; ++i) {
        if (matrix_argument_is_filtered(argc, argv, &i)) continue;
        child_args[child_argc++] = argv[i];
    }
    child_args[child_argc++] = physics_argument;
    child_args[child_argc++] = output_argument;
    child_args[child_argc] = NULL;
    int status = run_child(directory, result->binary, child_argc, child_args);
    free(child_args);
    return status;
}

static const char *matrix_status_label(int status) {
    if (status == 0) return "PASS";
    if (status == FPS_EXIT_GPU_CONTEXT_UNAVAILABLE ||
        status == FPS_EXIT_GPU_INITIALIZATION_FAILED) return "UNAVAILABLE";
    return "FAIL";
}

static bool write_matrix_report(const char *root, const char *scenario,
                                const MatrixResult *results, int result_count,
                                bool passed) {
    char path[MATRIX_PATH_CAPACITY];
    if (!join_path(path, sizeof(path), root, "matrix.json")) return false;
    FILE *file = fopen(path, "wb");
    if (!file) return false;
    fputs("{\n  \"schemaVersion\": 1,\n  \"scenario\": ", file);
    write_json_string(file, scenario);
    fprintf(file, ",\n  \"passed\": %s,\n  \"backends\": [\n", passed ? "true" : "false");
    for (int i = 0; i < result_count; ++i) {
        const MatrixResult *result = &results[i];
        fputs("    {\"requestedBackend\": ", file);
        write_json_string(file, result->backend);
        fputs(", \"status\": ", file);
        write_json_string(file, matrix_status_label(result->exit_status));
        fprintf(file, ", \"exitCode\": %d, \"report\": ", result->exit_status);
        write_json_string(file, result->report_path);
        fprintf(file, "}%s\n", (i + 1 < result_count) ? "," : "");
    }
    fputs("  ]\n}\n", file);
    int write_error = ferror(file);
    int close_error = fclose(file);
    if (!write_error && !close_error) {
        fprintf(stderr, "debug-matrix report=%s\n", path);
        return true;
    }
    return false;
}

static int run_debug_matrix(const char *directory, int argc, char **argv) {
    const char *scenario = option_value(argc, argv, "--debug-scenario=", "--debug-scenario");
    if (!scenario || !scenario[0]) {
        fprintf(stderr, "--debug-matrix requires --debug-scenario=<name>\n");
        return 2;
    }
    const char *requested_root = option_value(argc, argv, "--debug-output=", "--debug-output");
    char artifacts_root[MATRIX_PATH_CAPACITY];
    char matrix_root[MATRIX_PATH_CAPACITY];
    if (requested_root && requested_root[0]) {
        if (strlen(requested_root) >= sizeof(matrix_root)) return 5;
        strcpy(matrix_root, requested_root);
    } else {
        if (!join_path(artifacts_root, sizeof(artifacts_root), directory, "debug-artifacts") ||
            !join_path(matrix_root, sizeof(matrix_root), artifacts_root, scenario)) return 5;
    }
    if (!make_directory_recursive(matrix_root)) {
        fprintf(stderr, "Unable to create matrix output directory '%s': %s\n",
                matrix_root, strerror(errno));
        return 5;
    }

    MatrixResult results[] = {
        { .backend = "gpu", .binary = GPU_BINARY, .exit_status = -1 },
        { .backend = "cpu-mt", .binary = CPU_BINARY, .exit_status = -1 },
        { .backend = "cpu-st", .binary = CPU_BINARY, .exit_status = -1 }
    };
    int result_count = (int)(sizeof(results) / sizeof(results[0]));
    bool any_passed = false;
    bool any_failed = false;
    for (int i = 0; i < result_count; ++i) {
        const char *directory_name = (i == 0) ? "gpu-gl43" : results[i].backend;
        if (!join_path(results[i].output_dir, sizeof(results[i].output_dir),
                       matrix_root, directory_name) ||
            !join_path(results[i].report_path, sizeof(results[i].report_path),
                       results[i].output_dir, "report.json")) {
            results[i].exit_status = -1;
            any_failed = true;
            continue;
        }
        fprintf(stderr, "debug-matrix starting backend=%s output=%s\n",
                results[i].backend, results[i].output_dir);
        results[i].exit_status = run_matrix_child(directory, argc, argv, &results[i]);
        const char *label = matrix_status_label(results[i].exit_status);
        fprintf(stderr, "debug-matrix backend=%s status=%s exit=%d\n",
                results[i].backend, label, results[i].exit_status);
        if (results[i].exit_status == 0) any_passed = true;
        else if (strcmp(label, "UNAVAILABLE") != 0) any_failed = true;
    }
    bool passed = any_passed && !any_failed;
    if (!write_matrix_report(matrix_root, scenario, results, result_count, passed)) {
        fprintf(stderr, "Unable to write matrix report in '%s'\n", matrix_root);
        return 5;
    }
    return passed ? 0 : 5;
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
    if (matrix_requested(argc, argv)) {
        return run_debug_matrix(directory, argc, argv);
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
