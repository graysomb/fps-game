#!/usr/bin/env python3
"""
Interactive Call Graph Generator for fps_ray.c
Extracts all functions, docstrings, lines, caller-callee dependencies, and subsystems,
then generates an interactive, high-performance HTML visualization.
"""

import sys
import os
import re
import json

SUBSYSTEM_CONFIG = {
    "Game Loop & Lifecycle": {
        "color": "#ef4444",
        "badge_bg": "rgba(239, 68, 68, 0.18)",
        "badge_border": "#ef4444",
        "icon": "⚡",
        "description": "Engine startup, frame loops, main entrypoint, smoke tests and game resets"
    },
    "Physics / PBD": {
        "color": "#3b82f6",
        "badge_bg": "rgba(59, 130, 246, 0.18)",
        "badge_border": "#3b82f6",
        "icon": "⚛️",
        "description": "Position-Based Dynamics, particle integration, collision constraints, and VGS shape matching"
    },
    "Voxel Engine": {
        "color": "#10b981",
        "badge_bg": "rgba(16, 185, 129, 0.18)",
        "badge_border": "#10b981",
        "icon": "🧱",
        "description": "Spatial hash tables, voxel damage/carving, structural glue bonds, and surface caches"
    },
    "World Gen & Maps": {
        "color": "#84cc16",
        "badge_bg": "rgba(132, 204, 22, 0.18)",
        "badge_border": "#84cc16",
        "icon": "🗺️",
        "description": "Procedural arenas, fortresses, pillars, map save/load slots, and pickup spawning"
    },
    "Player & Combat": {
        "color": "#f59e0b",
        "badge_bg": "rgba(245, 158, 11, 0.18)",
        "badge_border": "#f59e0b",
        "icon": "🎯",
        "description": "Player motion, jumping, shooting, projectiles, melee swings, gravity tether, and matter damage"
    },
    "Bot AI": {
        "color": "#ec4899",
        "badge_bg": "rgba(236, 72, 153, 0.18)",
        "badge_border": "#ec4899",
        "icon": "🤖",
        "description": "Utility-based autonomous bot decision making (combat, harvest, flee, aim assist)"
    },
    "Networking": {
        "color": "#8b5cf6",
        "badge_bg": "rgba(139, 92, 246, 0.18)",
        "badge_border": "#8b5cf6",
        "icon": "🌐",
        "description": "LAN multiplayer, entity state serialization, client reconciliation, and voxel proxy synchronization"
    },
    "Rendering & UI": {
        "color": "#06b6d4",
        "badge_bg": "rgba(6, 182, 212, 0.18)",
        "badge_border": "#06b6d4",
        "icon": "🖥️",
        "description": "Split-screen camera layouts, greedy meshing, shaders, instanced cubes, HUD bars, and radar"
    },
    "Audio & SFX": {
        "color": "#f97316",
        "badge_bg": "rgba(249, 115, 22, 0.18)",
        "badge_border": "#f97316",
        "icon": "🔊",
        "description": "Sound synthesis, audio device lifecycle, and spatial sound playback"
    },
    "Creative Mode": {
        "color": "#14b8a6",
        "badge_bg": "rgba(20, 184, 166, 0.18)",
        "badge_border": "#14b8a6",
        "icon": "✏️",
        "description": "Noclip fly cameras, multi-voxel brush editing, and instant material palette placement"
    },
    "Math & Utilities": {
        "color": "#94a3b8",
        "badge_bg": "rgba(148, 163, 184, 0.18)",
        "badge_border": "#94a3b8",
        "icon": "📐",
        "description": "Vector3 math routines, clamping, color utilities, thread pool tasking, and bit helpers"
    }
}


def strip_comments_and_strings(src):
    chars = list(src)
    n = len(chars)
    i = 0
    while i < n:
        if chars[i] == "/" and i + 1 < n:
            if chars[i+1] == "/":
                chars[i] = " "
                chars[i+1] = " "
                i += 2
                while i < n and chars[i] != "\n":
                    chars[i] = " "
                    i += 1
                continue
            elif chars[i+1] == "*":
                chars[i] = " "
                chars[i+1] = " "
                i += 2
                while i < n - 1 and not (chars[i] == "*" and chars[i+1] == "/"):
                    if chars[i] != "\n":
                        chars[i] = " "
                    i += 1
                if i < n:
                    chars[i] = " "
                if i + 1 < n:
                    chars[i+1] = " "
                i += 2
                continue
        elif chars[i] == "\"":
            chars[i] = " "
            i += 1
            while i < n and chars[i] != "\"":
                if chars[i] == "\\":
                    chars[i] = " "
                    i += 1
                    if i < n and chars[i] != "\n":
                        chars[i] = " "
                elif chars[i] != "\n":
                    chars[i] = " "
                i += 1
            if i < n:
                chars[i] = " "
            i += 1
            continue
        elif chars[i] == "\'":
            chars[i] = " "
            i += 1
            while i < n and chars[i] != "\'":
                if chars[i] == "\\":
                    chars[i] = " "
                    i += 1
                    if i < n and chars[i] != "\n":
                        chars[i] = " "
                elif chars[i] != "\n":
                    chars[i] = " "
                i += 1
            if i < n:
                chars[i] = " "
            i += 1
            continue
        i += 1
    return "".join(chars)


def classify_subsystem(name, header, start_line):
    nl = name.lower()
    
    if "sfx" in nl or "sound" in nl or "audio" in nl:
        return "Audio & SFX"
        
    if "creative" in nl or ("brush" in nl and start_line > 13000):
        return "Creative Mode"
        
    if nl.startswith("net_") or "lan_" in nl or "packet" in nl or "proxy" in nl:
        return "Networking"
        
    if nl.startswith("bot_") or "bot" in nl or "calculateutility" in nl or nl in ("find_nearest_enemy", "init_bots", "update_bots"):
        return "Bot AI"
        
    if any(k in nl for k in ["buildstacked", "buildsplit", "buildleg", "buildbox", "buildfort", "carvefort", "buildgate", 
                             "buildtestworld", "buildprocedural", "buildblood", "builddebug", "builddemo", "pickup", "map_slot"]):
        return "World Gen & Maps"
        
    if any(k in nl for k in ["render", "draw", "hud", "mesh", "shader", "confetti", "instancing", "viewport", 
                             "view_layout", "visuals", "crosshair", "radar"]):
        return "Rendering & UI"
        
    if any(k in nl for k in ["player", "tether", "melee", "shoot", "bullet", "weapon", "jump", "damage", "matter",
                             "projectile", "smush", "keyboardinput", "gamepadinput", "aim_assist", "respawn", "kdratio", "perform_build"]):
        return "Player & Combat"
        
    if any(k in nl for k in ["pbd_", "physics_", "particle", "constraint", "vgs_", "ccd", "integrate", "wake_timer",
                             "break_mask", "pair_correction", "collision", "scratch", "belief", "cluster"]):
        return "Physics / PBD"
        
    if any(k in nl for k in ["voxel", "table_", "static_", "chunk", "glue", "hash", "face_visibility", "bounds", "grid"]):
        return "Voxel Engine"
        
    if (nl.startswith("v_") or nl.startswith("m_") or "clamp" in nl or "lerp" in nl or "vector" in nl or 
        "math" in nl or "bits" in nl or "cbrt" in nl or "pool" in nl or "thread" in nl or "ranges_overlap" in nl or "list_contains" in nl):
        return "Math & Utilities"
        
    if nl in ("main", "resetgame", "initgame", "rungame", "updategame", "drawgame", "cleanup") or "smoke" in nl:
        return "Game Loop & Lifecycle"
        
    if start_line < 2300:
        return "Physics / PBD"
    elif start_line < 3125:
        return "Math & Utilities"
    elif start_line < 6500:
        return "Voxel Engine"
    elif start_line < 7700:
        return "World Gen & Maps"
    elif start_line < 9300:
        return "Player & Combat"
    elif start_line < 12700:
        return "Physics / PBD"
    elif start_line < 15000:
        return "Rendering & UI"
    elif start_line < 15350:
        return "Player & Combat"
    elif start_line < 15650:
        return "Bot AI"
    elif start_line < 16210:
        return "Networking"
    else:
        return "Game Loop & Lifecycle"


def parse_functions(filepath):
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        source = f.read()

    source_lines = source.splitlines()
    cleaned = strip_comments_and_strings(source)
    
    functions = []
    paren_level = 0
    brace_level = 0
    header_start = 0

    for i, c in enumerate(cleaned):
        if c == "(":
            paren_level += 1
        elif c == ")":
            paren_level = max(0, paren_level - 1)
        elif c == "{" and paren_level == 0:
            if brace_level == 0:
                hdr = cleaned[header_start:i].strip()
                hdr_lines = [l.strip() for l in hdr.split("\n") if l.strip() and not l.strip().startswith("#")]
                hdr = " ".join(hdr_lines)
                if "(" in hdr and ")" in hdr and "=" not in hdr:
                    words = hdr.split()
                    if words and words[0] not in ("typedef", "struct", "enum", "union"):
                        m = re.search(r"([a-zA-Z_][a-zA-Z0-9_]*)\s*\([^()]*\)\s*$", hdr)
                        if m:
                            func_name = m.group(1)
                            if func_name not in ("if", "for", "while", "switch"):
                                start_line = source[:i].count("\n") + 1
                                functions.append({
                                    "name": func_name,
                                    "header": hdr,
                                    "start_idx": i,
                                    "start_line": start_line,
                                })
            brace_level += 1
        elif c == "}" and paren_level == 0:
            brace_level = max(0, brace_level - 1)
            if brace_level == 0:
                if functions and "end_line" not in functions[-1]:
                    functions[-1]["end_idx"] = i
                    functions[-1]["end_line"] = source[:i].count("\n") + 1
                header_start = i + 1
        elif c == ";" and brace_level == 0 and paren_level == 0:
            header_start = i + 1

    for f in functions:
        start_line = f["start_line"]
        end_line = f.get("end_line", start_line)
        f["line_count"] = end_line - start_line + 1
        
        doc_lines = []
        cur = start_line - 2
        while cur >= 0 and cur >= start_line - 7:
            line_str = source_lines[cur].strip()
            if line_str.startswith("//") or line_str.startswith("/*") or line_str.startswith("*"):
                cleaned_line = line_str.lstrip("/*").rstrip("*/").strip()
                if cleaned_line:
                    doc_lines.insert(0, cleaned_line)
                cur -= 1
            else:
                break
        f["doc"] = " ".join(doc_lines) if doc_lines else "Core implementation function in fps_ray.c"
        
        preview_lines = source_lines[start_line - 1 : min(len(source_lines), start_line + 9)]
        f["preview"] = "\n".join(preview_lines)
        f["subsystem"] = classify_subsystem(f["name"], f["header"], f["start_line"])

    return functions, cleaned, source


def extract_call_graph(functions, cleaned):
    func_map = {f["name"]: f for f in functions}
    func_names = set(func_map.keys())

    edges = []
    for f in functions:
        body = cleaned[f["start_idx"]:f["end_idx"]]
        calls = re.findall(r"\b([a-zA-Z_][a-zA-Z0-9_]*)\s*\(", body)
        internal_callees = set()
        external_calls = set()
        for c in calls:
            if c in func_names:
                if c != f["name"]:
                    internal_callees.add(c)
            elif not c.startswith("_") and c not in ("if", "for", "while", "switch", "sizeof", "return"):
                if any(c.startswith(p) for p in ("Draw", "rl", "Get", "Is", "Set", "Init", "Unload", "Vector3", "Matrix", "pthread", "malloc", "free", "memcpy", "memset")):
                    external_calls.add(c)

        f["callees"] = sorted(list(internal_callees))
        f["external_calls"] = sorted(list(external_calls))
        f["callers"] = []

    for f in functions:
        for callee_name in f["callees"]:
            if callee_name in func_map:
                func_map[callee_name]["callers"].append(f["name"])
                edges.append({
                    "source": f["name"],
                    "target": callee_name,
                    "cross_subsystem": f["subsystem"] != func_map[callee_name]["subsystem"]
                })

    for f in functions:
        f["callers"] = sorted(f["callers"])
        f["in_degree"] = len(f["callers"])
        f["out_degree"] = len(f["callees"])
        f["total_degree"] = f["in_degree"] + f["out_degree"]

    flow_matrix = {}
    for s1 in SUBSYSTEM_CONFIG:
        flow_matrix[s1] = {s2: 0 for s2 in SUBSYSTEM_CONFIG}

    for edge in edges:
        s_src = func_map[edge["source"]]["subsystem"]
        s_tgt = func_map[edge["target"]]["subsystem"]
        flow_matrix[s_src][s_tgt] += 1

    return edges, flow_matrix


def generate_html(functions, edges, flow_matrix, output_path):
    nodes_data = []
    for f in functions:
        nodes_data.append({
            "id": f["name"],
            "name": f["name"],
            "header": f["header"],
            "subsystem": f["subsystem"],
            "start_line": f["start_line"],
            "end_line": f["end_line"],
            "line_count": f["line_count"],
            "doc": f["doc"],
            "preview": f["preview"],
            "in_degree": f["in_degree"],
            "out_degree": f["out_degree"],
            "total_degree": f["total_degree"],
            "callers": f["callers"],
            "callees": f["callees"],
            "external_calls": f["external_calls"][:10],
        })

    subsystem_json = json.dumps(SUBSYSTEM_CONFIG)
    nodes_json = json.dumps(nodes_data)
    edges_json = json.dumps(edges)
    flow_json = json.dumps(flow_matrix)

    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>fps_ray.c Interactive Dependency Call Graph</title>
  <style>
    :root {
      --bg-dark: #0a0e17;
      --panel-bg: #111827;
      --card-bg: #1a2234;
      --border-color: #2e384d;
      --text-main: #f1f5f9;
      --text-muted: #94a3b8;
      --accent-cyan: #06b6d4;
      --accent-green: #10b981;
      --accent-blue: #3b82f6;
      --accent-amber: #f59e0b;
    }
    * {
      box-sizing: border-box;
      margin: 0;
      padding: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      user-select: none;
    }
    body {
      background-color: var(--bg-dark);
      color: var(--text-main);
      overflow: hidden;
      height: 100vh;
      display: flex;
      flex-direction: column;
    }
    header {
      height: 60px;
      background-color: var(--panel-bg);
      border-bottom: 1px solid var(--border-color);
      display: flex;
      align-items: center;
      justify-content: space-between;
      padding: 0 20px;
      z-index: 10;
    }
    .brand {
      display: flex;
      align-items: center;
      gap: 12px;
    }
    .brand-icon {
      width: 32px;
      height: 32px;
      background: linear-gradient(135deg, #06b6d4, #3b82f6);
      border-radius: 8px;
      display: flex;
      align-items: center;
      justify-content: center;
      font-weight: bold;
      font-size: 18px;
      color: #fff;
    }
    .brand-text h1 {
      font-size: 16px;
      font-weight: 700;
      letter-spacing: 0.5px;
      color: #fff;
    }
    .brand-text p {
      font-size: 11px;
      color: var(--text-muted);
    }
    .header-actions {
      display: flex;
      align-items: center;
      gap: 12px;
    }
    .search-wrapper {
      position: relative;
      width: 280px;
    }
    .search-input {
      width: 100%;
      background: #0a0e17;
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 8px 12px 8px 34px;
      color: #fff;
      font-size: 13px;
      outline: none;
      transition: border-color 0.2s;
    }
    .search-input:focus {
      border-color: var(--accent-cyan);
    }
    .search-icon {
      position: absolute;
      left: 10px;
      top: 50%;
      transform: translateY(-50%);
      color: var(--text-muted);
      font-size: 14px;
    }
    .autocomplete-list {
      position: absolute;
      top: 100%;
      left: 0;
      right: 0;
      background: var(--panel-bg);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      margin-top: 4px;
      max-height: 280px;
      overflow-y: auto;
      display: none;
      z-index: 50;
      box-shadow: 0 8px 24px rgba(0,0,0,0.5);
    }
    .autocomplete-item {
      padding: 8px 12px;
      font-size: 12px;
      cursor: pointer;
      display: flex;
      align-items: center;
      justify-content: space-between;
      border-bottom: 1px solid rgba(255,255,255,0.05);
    }
    .autocomplete-item:hover {
      background: #1e293b;
    }
    .badge {
      font-size: 10px;
      padding: 2px 6px;
      border-radius: 4px;
      font-weight: 600;
    }
    .btn {
      background: #1e293b;
      border: 1px solid var(--border-color);
      color: var(--text-main);
      padding: 6px 14px;
      border-radius: 6px;
      font-size: 12px;
      font-weight: 500;
      cursor: pointer;
      display: flex;
      align-items: center;
      gap: 6px;
      transition: all 0.15s;
    }
    .btn:hover {
      background: #2b3952;
      border-color: var(--accent-cyan);
    }
    .btn.active {
      background: var(--accent-cyan);
      color: #000;
      font-weight: 600;
      border-color: var(--accent-cyan);
    }
    .main-content {
      flex: 1;
      display: flex;
      position: relative;
      overflow: hidden;
    }
    #graph-container {
      flex: 1;
      position: relative;
      overflow: hidden;
      background: radial-gradient(circle at center, #111827 0%, #06090e 100%);
    }
    canvas {
      display: block;
      width: 100%;
      height: 100%;
      cursor: grab;
    }
    canvas:active {
      cursor: grabbing;
    }
    .floating-hud {
      position: absolute;
      top: 16px;
      left: 16px;
      display: flex;
      flex-direction: column;
      gap: 8px;
      z-index: 5;
    }
    .preset-bar {
      display: flex;
      align-items: center;
      background: rgba(17, 24, 39, 0.85);
      backdrop-filter: blur(8px);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      padding: 4px;
      gap: 4px;
    }
    .degree-slider-container {
      display: flex;
      align-items: center;
      gap: 8px;
      margin-left: 12px;
      padding-left: 12px;
      border-left: 1px solid var(--border-color);
      font-size: 11px;
      color: var(--text-muted);
    }
    .degree-slider {
      accent-color: var(--accent-cyan);
      cursor: pointer;
      width: 90px;
    }
    .subsystem-filter-bar {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      max-width: 820px;
      background: rgba(17, 24, 39, 0.85);
      backdrop-filter: blur(8px);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      padding: 8px;
    }
    .sub-pill {
      font-size: 11px;
      padding: 4px 8px;
      border-radius: 12px;
      cursor: pointer;
      display: flex;
      align-items: center;
      gap: 6px;
      border: 1px solid transparent;
      opacity: 0.45;
      transition: all 0.15s;
    }
    .sub-pill.active {
      opacity: 1;
    }
    .sub-pill-dot {
      width: 8px;
      height: 8px;
      border-radius: 50%;
    }
    .zoom-controls {
      position: absolute;
      bottom: 24px;
      left: 20px;
      display: flex;
      flex-direction: column;
      background: rgba(17, 24, 39, 0.85);
      backdrop-filter: blur(8px);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      overflow: hidden;
      z-index: 5;
    }
    .zoom-btn {
      width: 36px;
      height: 36px;
      display: flex;
      align-items: center;
      justify-content: center;
      background: transparent;
      border: none;
      color: #fff;
      font-size: 16px;
      cursor: pointer;
      border-bottom: 1px solid var(--border-color);
    }
    .zoom-btn:last-child {
      border-bottom: none;
    }
    .zoom-btn:hover {
      background: #1e293b;
      color: var(--accent-cyan);
    }
    .legend-card {
      position: absolute;
      bottom: 24px;
      right: 440px;
      background: rgba(17, 24, 39, 0.85);
      backdrop-filter: blur(8px);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      padding: 10px 14px;
      font-size: 11px;
      z-index: 5;
      display: flex;
      gap: 16px;
      align-items: center;
    }
    .legend-item {
      display: flex;
      align-items: center;
      gap: 6px;
    }
    .legend-arrow {
      width: 20px;
      height: 2px;
      position: relative;
    }
    .legend-arrow.caller {
      background: var(--accent-green);
    }
    .legend-arrow.callee {
      background: var(--accent-cyan);
    }
    .inspector-panel {
      width: 420px;
      background: var(--panel-bg);
      border-left: 1px solid var(--border-color);
      display: flex;
      flex-direction: column;
      height: 100%;
      z-index: 10;
      box-shadow: -4px 0 24px rgba(0,0,0,0.4);
      transition: transform 0.25s ease;
    }
    .inspector-header {
      padding: 16px 20px;
      border-bottom: 1px solid var(--border-color);
      background: #131d2e;
    }
    .inspector-header .badge-row {
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 8px;
    }
    .func-name {
      font-size: 18px;
      font-weight: 700;
      color: #fff;
      word-break: break-all;
      font-family: monospace;
    }
    .func-loc {
      font-size: 12px;
      color: var(--accent-cyan);
      margin-top: 4px;
      display: flex;
      align-items: center;
      gap: 6px;
    }
    .func-loc a {
      color: inherit;
      text-decoration: none;
    }
    .func-loc a:hover {
      text-decoration: underline;
    }
    .stats-row {
      display: grid;
      grid-template-columns: repeat(3, 1fr);
      gap: 8px;
      margin-top: 14px;
    }
    .stat-box {
      background: var(--card-bg);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 8px;
      text-align: center;
    }
    .stat-val {
      font-size: 16px;
      font-weight: 700;
      color: #fff;
    }
    .stat-label {
      font-size: 10px;
      color: var(--text-muted);
      text-transform: uppercase;
      margin-top: 2px;
    }
    .inspector-body {
      flex: 1;
      overflow-y: auto;
      padding: 16px 20px;
      display: flex;
      flex-direction: column;
      gap: 16px;
    }
    .section-title {
      font-size: 11px;
      text-transform: uppercase;
      letter-spacing: 0.8px;
      color: var(--text-muted);
      font-weight: 700;
      margin-bottom: 8px;
      display: flex;
      justify-content: space-between;
      align-items: center;
    }
    .doc-box {
      background: var(--card-bg);
      border-left: 3px solid var(--accent-blue);
      border-radius: 4px;
      padding: 10px 12px;
      font-size: 12px;
      color: #cbd5e1;
      line-height: 1.5;
    }
    .code-preview {
      background: #080c14;
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 10px;
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
      font-size: 11px;
      color: #e2e8f0;
      overflow-x: auto;
      white-space: pre;
      line-height: 1.4;
      max-height: 180px;
    }
    .links-list {
      display: flex;
      flex-direction: column;
      gap: 6px;
      max-height: 160px;
      overflow-y: auto;
    }
    .link-item {
      background: var(--card-bg);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 6px 10px;
      font-size: 12px;
      font-family: monospace;
      cursor: pointer;
      display: flex;
      align-items: center;
      justify-content: space-between;
      transition: all 0.15s;
    }
    .link-item:hover {
      background: #243048;
      border-color: var(--accent-cyan);
      transform: translateX(3px);
    }
    .empty-state {
      font-size: 12px;
      color: var(--text-muted);
      font-style: italic;
      padding: 8px 0;
    }
    #tooltip {
      position: absolute;
      background: rgba(15, 23, 42, 0.95);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 8px 12px;
      font-size: 12px;
      pointer-events: none;
      z-index: 100;
      display: none;
      box-shadow: 0 4px 16px rgba(0,0,0,0.4);
    }
    #tooltip .tt-title {
      font-weight: 700;
      font-family: monospace;
      color: #fff;
      margin-bottom: 4px;
    }
    #tooltip .tt-sub {
      font-size: 11px;
      color: var(--text-muted);
    }
    #flow-modal {
      position: absolute;
      top: 60px;
      left: 0;
      right: 0;
      bottom: 0;
      background: rgba(10, 14, 23, 0.92);
      backdrop-filter: blur(12px);
      z-index: 40;
      display: none;
      padding: 30px;
      overflow-y: auto;
    }
    .flow-modal-inner {
      max-width: 1000px;
      margin: 0 auto;
      background: var(--panel-bg);
      border: 1px solid var(--border-color);
      border-radius: 12px;
      padding: 24px;
    }
    .modal-header {
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 20px;
    }
    .modal-title {
      font-size: 18px;
      font-weight: 700;
    }
    .flow-table {
      width: 100%;
      border-collapse: collapse;
      font-size: 11px;
      margin-top: 14px;
    }
    .flow-table th, .flow-table td {
      padding: 8px 10px;
      border: 1px solid var(--border-color);
      text-align: center;
    }
    .flow-table th {
      background: #1a2234;
      font-weight: 600;
      white-space: nowrap;
    }
    .flow-table td.flow-cell {
      background: #0d131f;
      font-weight: 600;
    }
    .flow-table td.flow-cell.has-calls {
      background: rgba(6, 182, 212, 0.15);
      color: var(--accent-cyan);
    }
    .close-btn {
      background: none;
      border: none;
      color: #fff;
      font-size: 20px;
      cursor: pointer;
    }
  </style>
</head>
<body>

  <header>
    <div class="brand">
      <div class="brand-icon">⚡</div>
      <div class="brand-text">
        <h1>fps_ray.c Function Dependency Graph</h1>
        <p>17,310 Lines &bull; """ + str(len(functions)) + """ Functions &bull; """ + str(len(edges)) + """ Call Dependencies</p>
      </div>
    </div>

    <div class="header-actions">
      <div class="search-wrapper">
        <span class="search-icon">🔍</span>
        <input type="text" id="search-input" class="search-input" placeholder="Search function (or press /)...">
        <div id="autocomplete" class="autocomplete-list"></div>
      </div>

      <button id="btn-flow" class="btn" title="View subsystem dependency matrix">📊 Subsystem Flow</button>
      <button id="btn-labels" class="btn" title="Toggle label density">🏷️ All Labels</button>
      <button id="btn-physics-toggle" class="btn" title="Freeze or resume physics forces">⏸️ Freeze</button>
      <button id="btn-reset" class="btn">↺ Reset View</button>
    </div>
  </header>

  <div class="main-content">
    <div id="graph-container">
      <canvas id="graph-canvas"></canvas>

      <div class="floating-hud">
        <div class="preset-bar">
          <button class="btn active" data-preset="all">Full (""" + str(len(functions)) + """)</button>
          <button class="btn" data-preset="core">Core Hubs</button>
          <button class="btn" data-preset="physics">Physics Pipeline</button>
          <button class="btn" data-preset="voxel">Voxel Engine</button>
          <button class="btn" data-preset="net">Net & Sync</button>
          <button class="btn" data-preset="combat">Combat</button>

          <div class="degree-slider-container">
            <span>Min Calls:</span>
            <input type="range" id="degree-slider" class="degree-slider" min="0" max="25" value="0">
            <span id="degree-val" style="color:#fff;font-weight:700;">0</span>
            <span id="visible-count-badge" style="color:var(--accent-cyan);">(""" + str(len(functions)) + """)</span>
          </div>
        </div>

        <div id="subsystem-filters" class="subsystem-filter-bar"></div>
      </div>

      <div class="zoom-controls">
        <button id="zoom-in" class="zoom-btn" title="Zoom In">+</button>
        <button id="zoom-out" class="zoom-btn" title="Zoom Out">&minus;</button>
        <button id="zoom-fit" class="zoom-btn" title="Fit All Nodes">⛶</button>
      </div>

      <div class="legend-card">
        <div class="legend-item">
          <div class="legend-arrow caller"></div>
          <span>Caller (Incoming)</span>
        </div>
        <div class="legend-item">
          <div class="legend-arrow callee"></div>
          <span>Callee (Outgoing)</span>
        </div>
        <div class="legend-item">
          <span style="color: var(--text-muted);">&bull; Node size = Total call degree</span>
        </div>
      </div>
    </div>

    <aside id="inspector" class="inspector-panel">
      <div class="inspector-header">
        <div class="badge-row">
          <span id="insp-subsystem-badge" class="badge">Subsystem</span>
          <span id="insp-scope" style="font-size: 11px; color: var(--text-muted);">Source Function</span>
        </div>
        <div id="insp-name" class="func-name">function_name</div>
        <div class="func-loc">
          <span>📍</span>
          <a id="insp-line-link" href="#" target="_blank">fps_ray.c: L0 - L0</a>
        </div>

        <div class="stats-row">
          <div class="stat-box">
            <div id="insp-callers-count" class="stat-val" style="color: var(--accent-green);">0</div>
            <div class="stat-label">Callers</div>
          </div>
          <div class="stat-box">
            <div id="insp-callees-count" class="stat-val" style="color: var(--accent-cyan);">0</div>
            <div class="stat-label">Callees</div>
          </div>
          <div class="stat-box">
            <div id="insp-lines-count" class="stat-val" style="color: var(--accent-amber);">0</div>
            <div class="stat-label">Lines of C</div>
          </div>
        </div>
      </div>

      <div class="inspector-body">
        <div>
          <div class="section-title">Summary / Documentation</div>
          <div id="insp-doc" class="doc-box">Function documentation.</div>
        </div>

        <div>
          <div class="section-title">Code Preview</div>
          <pre id="insp-code" class="code-preview"></pre>
        </div>

        <div>
          <div class="section-title">
            <span>Callers (Called By)</span>
            <span id="insp-callers-badge" style="color: var(--accent-green);">0</span>
          </div>
          <div id="insp-callers-list" class="links-list"></div>
        </div>

        <div>
          <div class="section-title">
            <span>Callees (Calls)</span>
            <span id="insp-callees-badge" style="color: var(--accent-cyan);">0</span>
          </div>
          <div id="insp-callees-list" class="links-list"></div>
        </div>

        <div>
          <div class="section-title">
            <span>External Library Calls</span>
            <span id="insp-external-badge" style="color: var(--text-muted);">0</span>
          </div>
          <div id="insp-external-list" class="links-list"></div>
        </div>
      </div>
    </aside>
  </div>

  <div id="flow-modal">
    <div class="flow-modal-inner">
      <div class="modal-header">
        <div>
          <div class="modal-title">Subsystem Dependency Matrix</div>
          <p style="font-size: 12px; color: var(--text-muted); margin-top: 4px;">Shows direct call volumes between architectural subsystems (Row = Caller Subsystem, Column = Callee Subsystem)</p>
        </div>
        <button id="close-flow-modal" class="close-btn">&times;</button>
      </div>
      <div style="overflow-x: auto;">
        <table id="flow-table-content" class="flow-table"></table>
      </div>
    </div>
  </div>

  <div id="tooltip">
    <div id="tt-name" class="tt-title"></div>
    <div id="tt-sub" class="tt-sub"></div>
    <div id="tt-stats" style="font-size: 11px; margin-top: 4px; color: #cbd5e1;"></div>
  </div>

  <script>
    const SUBSYSTEMS = """ + subsystem_json + """;
    const RAW_NODES = """ + nodes_json + """;
    const RAW_EDGES = """ + edges_json + """;
    const FLOW_MATRIX = """ + flow_json + """;

    const state = {
      nodes: [],
      edges: [],
      nodeMap: new Map(),
      activeSubsystems: new Set(Object.keys(SUBSYSTEMS)),
      selectedNode: null,
      hoveredNode: null,
      zoom: 0.85,
      panX: 0,
      panY: 0,
      isDragging: false,
      draggedNode: null,
      lastMouseX: 0,
      lastMouseY: 0,
      physicsRunning: true,
      currentPreset: "all",
      minDegreeFilter: 0,
      showAllLabels: false
    };

    const canvas = document.getElementById("graph-canvas");
    const ctx = canvas.getContext("2d");
    let width, height;

    function resizeCanvas() {
      const container = document.getElementById("graph-container");
      width = container.clientWidth;
      height = container.clientHeight;
      const dpr = window.devicePixelRatio || 1;
      canvas.width = width * dpr;
      canvas.height = height * dpr;
      ctx.scale(dpr, dpr);
    }
    window.addEventListener("resize", () => {
      resizeCanvas();
      draw();
    });

    function initNodesAndPositions() {
      const subKeys = Object.keys(SUBSYSTEMS);
      const subAngles = {};
      subKeys.forEach((key, idx) => {
        subAngles[key] = (idx / subKeys.length) * Math.PI * 2;
      });

      state.nodes = RAW_NODES.map((n) => {
        const angle = subAngles[n.subsystem] || 0;
        const radius = 350 + Math.random() * 400;
        const jitterX = (Math.random() - 0.5) * 200;
        const jitterY = (Math.random() - 0.5) * 200;

        const node = {
          ...n,
          x: Math.cos(angle) * radius + jitterX,
          y: Math.sin(angle) * radius + jitterY,
          vx: 0,
          vy: 0,
          radius: Math.max(5, Math.min(22, 5 + Math.sqrt(n.total_degree) * 2.8)),
          color: SUBSYSTEMS[n.subsystem]?.color || "#94a3b8",
          visible: true
        };
        state.nodeMap.set(n.id, node);
        return node;
      });

      state.edges = RAW_EDGES.map(e => ({
        source: state.nodeMap.get(e.source),
        target: state.nodeMap.get(e.target),
        cross_subsystem: e.cross_subsystem
      })).filter(e => e.source && e.target);

      state.panX = width / 2;
      state.panY = height / 2;
    }

    function stepPhysics() {
      if (!state.physicsRunning) return;

      const nodes = state.nodes.filter(n => n.visible);
      const kRepulsion = 1400;
      const kSpring = 0.0035;
      const springLength = 80;
      const kCentering = 0.0008;

      for (let i = 0; i < nodes.length; i++) {
        const n = nodes[i];
        if (n === state.draggedNode) continue;
        n.vx -= n.x * kCentering;
        n.vy -= n.y * kCentering;
      }

      for (let i = 0; i < nodes.length; i++) {
        const n1 = nodes[i];
        for (let j = i + 1; j < nodes.length; j++) {
          const n2 = nodes[j];
          const dx = n2.x - n1.x;
          const dy = n2.y - n1.y;
          const distSq = dx * dx + dy * dy + 10;
          if (distSq < 160000) {
            const dist = Math.sqrt(distSq);
            const force = kRepulsion / distSq;
            const fx = (dx / dist) * force;
            const fy = (dy / dist) * force;

            if (n1 !== state.draggedNode) { n1.vx -= fx; n1.vy -= fy; }
            if (n2 !== state.draggedNode) { n2.vx += fx; n2.vy += fy; }
          }
        }
      }

      for (let i = 0; i < state.edges.length; i++) {
        const e = state.edges[i];
        if (!e.source.visible || !e.target.visible) continue;

        const dx = e.target.x - e.source.x;
        const dy = e.target.y - e.source.y;
        const dist = Math.sqrt(dx * dx + dy * dy) || 1;
        const force = (dist - springLength) * kSpring;
        const fx = (dx / dist) * force;
        const fy = (dy / dist) * force;

        if (e.source !== state.draggedNode) { e.source.vx += fx; e.source.vy += fy; }
        if (e.target !== state.draggedNode) { e.target.vx -= fx; e.target.vy -= fy; }
      }

      const damping = 0.88;
      for (let i = 0; i < nodes.length; i++) {
        const n = nodes[i];
        if (n === state.draggedNode) continue;
        n.vx *= damping;
        n.vy *= damping;
        n.x += n.vx;
        n.y += n.vy;
      }
    }

    function draw() {
      stepPhysics();

      ctx.clearRect(0, 0, width, height);

      ctx.save();
      ctx.translate(state.panX, state.panY);
      ctx.scale(state.zoom, state.zoom);

      const sel = state.selectedNode;
      const hov = state.hoveredNode;
      const activeFocus = sel || hov;

      const directCallers = activeFocus ? new Set(activeFocus.callers) : null;
      const directCallees = activeFocus ? new Set(activeFocus.callees) : null;

      for (let i = 0; i < state.edges.length; i++) {
        const e = state.edges[i];
        if (!e.source.visible || !e.target.visible) continue;

        let strokeColor = "rgba(255, 255, 255, 0.05)";
        let lineWidth = 1;
        let isHighlighted = false;

        if (activeFocus) {
          if (e.target.id === activeFocus.id && directCallers.has(e.source.id)) {
            strokeColor = "#10b981";
            lineWidth = 2.2;
            isHighlighted = true;
          } else if (e.source.id === activeFocus.id && directCallees.has(e.target.id)) {
            strokeColor = "#06b6d4";
            lineWidth = 2.2;
            isHighlighted = true;
          } else {
            strokeColor = "rgba(255, 255, 255, 0.02)";
          }
        }

        ctx.strokeStyle = strokeColor;
        ctx.lineWidth = lineWidth;
        ctx.beginPath();
        ctx.moveTo(e.source.x, e.source.y);
        ctx.lineTo(e.target.x, e.target.y);
        ctx.stroke();

        if (isHighlighted) {
          const dx = e.target.x - e.source.x;
          const dy = e.target.y - e.source.y;
          const len = Math.sqrt(dx * dx + dy * dy);
          if (len > 20) {
            const mx = (e.source.x + e.target.x) / 2;
            const my = (e.source.y + e.target.y) / 2;
            const angle = Math.atan2(dy, dx);
            const arrowSize = 7;

            ctx.fillStyle = strokeColor;
            ctx.beginPath();
            ctx.moveTo(mx + Math.cos(angle) * arrowSize, my + Math.sin(angle) * arrowSize);
            ctx.lineTo(mx + Math.cos(angle + 2.5) * arrowSize, my + Math.sin(angle + 2.5) * arrowSize);
            ctx.lineTo(mx + Math.cos(angle - 2.5) * arrowSize, my + Math.sin(angle - 2.5) * arrowSize);
            ctx.closePath();
            ctx.fill();
          }
        }
      }

      for (let i = 0; i < state.nodes.length; i++) {
        const n = state.nodes[i];
        if (!n.visible) continue;

        let alpha = 1.0;
        let isFocused = false;
        let isCaller = false;
        let isCallee = false;

        if (activeFocus) {
          if (n.id === activeFocus.id) {
            isFocused = true;
          } else if (directCallers.has(n.id)) {
            isCaller = true;
          } else if (directCallees.has(n.id)) {
            isCallee = true;
          } else {
            alpha = 0.12;
          }
        }

        ctx.save();
        ctx.globalAlpha = alpha;

        ctx.beginPath();
        ctx.arc(n.x, n.y, n.radius, 0, Math.PI * 2);

        if (isFocused) {
          ctx.fillStyle = "#ffffff";
          ctx.fill();
          ctx.lineWidth = 4;
          ctx.strokeStyle = "#06b6d4";
          ctx.stroke();

          ctx.beginPath();
          ctx.arc(n.x, n.y, n.radius + 8, 0, Math.PI * 2);
          ctx.strokeStyle = "rgba(6, 182, 212, 0.45)";
          ctx.lineWidth = 2;
          ctx.stroke();
        } else if (isCaller) {
          ctx.fillStyle = "#10b981";
          ctx.fill();
          ctx.lineWidth = 3;
          ctx.strokeStyle = "#ffffff";
          ctx.stroke();
        } else if (isCallee) {
          ctx.fillStyle = "#06b6d4";
          ctx.fill();
          ctx.lineWidth = 3;
          ctx.strokeStyle = "#ffffff";
          ctx.stroke();
        } else {
          ctx.fillStyle = n.color;
          ctx.fill();
          ctx.lineWidth = 1.5;
          ctx.strokeStyle = "rgba(255,255,255,0.3)";
          ctx.stroke();
        }

        const shouldShowLabel = state.showAllLabels || isFocused || isCaller || isCallee || (state.zoom > 1.4) || (n.total_degree >= 18);
        if (shouldShowLabel) {
          ctx.font = isFocused ? "bold 13px monospace" : "11px monospace";
          ctx.textAlign = "center";
          ctx.textBaseline = "middle";

          ctx.strokeStyle = "#0a0e17";
          ctx.lineWidth = 3;
          ctx.strokeText(n.name, n.x, n.y + n.radius + 12);

          ctx.fillStyle = isFocused ? "#ffffff" : "#cbd5e1";
          ctx.fillText(n.name, n.x, n.y + n.radius + 12);
        }

        ctx.restore();
      }

      ctx.restore();
      requestAnimationFrame(draw);
    }

    function screenToWorld(sx, sy) {
      return {
        x: (sx - state.panX) / state.zoom,
        y: (sy - state.panY) / state.zoom
      };
    }

    function findNodeAt(worldX, worldY) {
      for (let i = state.nodes.length - 1; i >= 0; i--) {
        const n = state.nodes[i];
        if (!n.visible) continue;
        const dx = n.x - worldX;
        const dy = n.y - worldY;
        if (dx * dx + dy * dy <= (n.radius + 4) * (n.radius + 4)) {
          return n;
        }
      }
      return null;
    }

    canvas.addEventListener("mousedown", (e) => {
      const rect = canvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;
      const w = screenToWorld(mx, my);

      const clickedNode = findNodeAt(w.x, w.y);
      if (clickedNode) {
        state.draggedNode = clickedNode;
        selectNode(clickedNode);
      } else {
        state.isDragging = true;
        state.lastMouseX = e.clientX;
        state.lastMouseY = e.clientY;
      }
    });

    window.addEventListener("mousemove", (e) => {
      const rect = canvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;

      if (state.draggedNode) {
        const w = screenToWorld(mx, my);
        state.draggedNode.x = w.x;
        state.draggedNode.y = w.y;
        state.draggedNode.vx = 0;
        state.draggedNode.vy = 0;
      } else if (state.isDragging) {
        state.panX += e.clientX - state.lastMouseX;
        state.panY += e.clientY - state.lastMouseY;
        state.lastMouseX = e.clientX;
        state.lastMouseY = e.clientY;
      } else {
        const w = screenToWorld(mx, my);
        const hovered = findNodeAt(w.x, w.y);
        state.hoveredNode = hovered;

        const tt = document.getElementById("tooltip");
        if (hovered) {
          tt.style.display = "block";
          tt.style.left = (e.clientX + 14) + "px";
          tt.style.top = (e.clientY + 14) + "px";
          document.getElementById("tt-name").textContent = hovered.name + "()";
          document.getElementById("tt-sub").textContent = hovered.subsystem + " • L" + hovered.start_line + "-L" + hovered.end_line;
          document.getElementById("tt-stats").textContent = `Callers: ${hovered.in_degree} | Callees: ${hovered.out_degree} | LOC: ${hovered.line_count}`;
        } else {
          tt.style.display = "none";
        }
      }
    });

    window.addEventListener("mouseup", () => {
      state.isDragging = false;
      state.draggedNode = null;
    });

    canvas.addEventListener("wheel", (e) => {
      e.preventDefault();
      const zoomFactor = e.deltaY < 0 ? 1.12 : 0.88;
      const rect = canvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;

      const wBefore = screenToWorld(mx, my);
      state.zoom = Math.max(0.15, Math.min(4.0, state.zoom * zoomFactor));
      const wAfter = screenToWorld(mx, my);

      state.panX += (wAfter.x - wBefore.x) * state.zoom;
      state.panY += (wAfter.y - wBefore.y) * state.zoom;
    }, { passive: false });

    canvas.addEventListener("dblclick", (e) => {
      const rect = canvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;
      const w = screenToWorld(mx, my);
      const clicked = findNodeAt(w.x, w.y);
      if (!clicked) {
        clearSelection();
      }
    });

    function selectNode(node) {
      state.selectedNode = node;
      if (!node) return;

      document.getElementById("insp-name").textContent = node.name + "()";
      const badge = document.getElementById("insp-subsystem-badge");
      badge.textContent = node.subsystem;
      badge.style.backgroundColor = SUBSYSTEMS[node.subsystem]?.badge_bg || "#334155";
      badge.style.color = SUBSYSTEMS[node.subsystem]?.color || "#fff";
      badge.style.border = "1px solid " + (SUBSYSTEMS[node.subsystem]?.badge_border || "transparent");

      const lineLink = document.getElementById("insp-line-link");
      lineLink.textContent = `fps_ray.c: L${node.start_line} - L${node.end_line} (${node.line_count} lines)`;
      lineLink.href = `fps_ray.c#L${node.start_line}`;

      document.getElementById("insp-callers-count").textContent = node.in_degree;
      document.getElementById("insp-callees-count").textContent = node.out_degree;
      document.getElementById("insp-lines-count").textContent = node.line_count;

      document.getElementById("insp-doc").textContent = node.doc;
      document.getElementById("insp-code").textContent = node.preview;

      const callersList = document.getElementById("insp-callers-list");
      callersList.innerHTML = "";
      document.getElementById("insp-callers-badge").textContent = node.callers.length;
      if (node.callers.length === 0) {
        callersList.innerHTML = '<div class="empty-state">Root / Entry function (no internal callers)</div>';
      } else {
        node.callers.forEach(cName => {
          const item = document.createElement("div");
          item.className = "link-item";
          item.innerHTML = `<span>⚡ ${cName}()</span><span style="font-size: 10px; color: var(--accent-green);">CALLER</span>`;
          item.onclick = () => jumpToNode(cName);
          callersList.appendChild(item);
        });
      }

      const calleesList = document.getElementById("insp-callees-list");
      calleesList.innerHTML = "";
      document.getElementById("insp-callees-badge").textContent = node.callees.length;
      if (node.callees.length === 0) {
        callersList.innerHTML = '<div class="empty-state">Leaf function (no internal outgoing calls)</div>';
      } else {
        node.callees.forEach(cName => {
          const item = document.createElement("div");
          item.className = "link-item";
          item.innerHTML = `<span>⚡ ${cName}()</span><span style="font-size: 10px; color: var(--accent-cyan);">CALLEE</span>`;
          item.onclick = () => jumpToNode(cName);
          calleesList.appendChild(item);
        });
      }

      const extList = document.getElementById("insp-external-list");
      extList.innerHTML = "";
      document.getElementById("insp-external-badge").textContent = node.external_calls.length;
      if (node.external_calls.length === 0) {
        extList.innerHTML = '<div class="empty-state">No external API calls</div>';
      } else {
        node.external_calls.forEach(eName => {
          const item = document.createElement("div");
          item.className = "link-item";
          item.innerHTML = `<span>📦 ${eName}()</span><span style="font-size: 10px; color: var(--text-muted);">LIB</span>`;
          extList.appendChild(item);
        });
      }
    }

    function clearSelection() {
      state.selectedNode = null;
    }

    function jumpToNode(nodeId) {
      const node = state.nodeMap.get(nodeId);
      if (!node) return;

      if (!state.activeSubsystems.has(node.subsystem)) {
        state.activeSubsystems.add(node.subsystem);
        applyFilters();
      }

      selectNode(node);

      state.panX = width / 2 - node.x * state.zoom;
      state.panY = height / 2 - node.y * state.zoom;
    }

    function setupSubsystemFilters() {
      const filterContainer = document.getElementById("subsystem-filters");
      filterContainer.innerHTML = "";

      Object.entries(SUBSYSTEMS).forEach(([name, cfg]) => {
        const count = RAW_NODES.filter(n => n.subsystem === name).length;
        const pill = document.createElement("div");
        pill.className = "sub-pill active";
        pill.dataset.sub = name;
        pill.style.backgroundColor = cfg.badge_bg;
        pill.style.borderColor = cfg.badge_border;
        pill.style.color = "#fff";
        pill.innerHTML = `<span class="sub-pill-dot" style="background:${cfg.color};"></span><span>${cfg.icon} ${name} (${count})</span>`;

        pill.onclick = (e) => {
          if (e.altKey) {
            state.activeSubsystems.clear();
            state.activeSubsystems.add(name);
          } else {
            if (state.activeSubsystems.has(name)) {
              if (state.activeSubsystems.size > 1) state.activeSubsystems.delete(name);
            } else {
              state.activeSubsystems.add(name);
            }
          }
          applyFilters();
        };
        filterContainer.appendChild(pill);
      });
    }

    function applyFilters() {
      document.querySelectorAll(".sub-pill").forEach(pill => {
        const sub = pill.dataset.sub;
        if (state.activeSubsystems.has(sub)) {
          pill.classList.add("active");
        } else {
          pill.classList.remove("active");
        }
      });

      let visibleCount = 0;
      state.nodes.forEach(n => {
        let visible = state.activeSubsystems.has(n.subsystem);
        if (n.total_degree < state.minDegreeFilter) visible = false;
        if (state.currentPreset === "core" && n.total_degree < 10) visible = false;
        n.visible = visible;
        if (visible) visibleCount++;
      });

      document.getElementById("visible-count-badge").textContent = `(${visibleCount})`;
    }

    document.querySelectorAll(".preset-bar .btn").forEach(btn => {
      btn.onclick = () => {
        document.querySelectorAll(".preset-bar .btn").forEach(b => b.classList.remove("active"));
        btn.classList.add("active");
        const preset = btn.dataset.preset;
        state.currentPreset = preset;

        if (preset === "all") {
          state.activeSubsystems = new Set(Object.keys(SUBSYSTEMS));
          state.minDegreeFilter = 0;
        } else if (preset === "core") {
          state.activeSubsystems = new Set(Object.keys(SUBSYSTEMS));
          state.minDegreeFilter = 10;
        } else if (preset === "physics") {
          state.activeSubsystems = new Set(["Physics / PBD", "Math & Utilities"]);
          state.minDegreeFilter = 0;
        } else if (preset === "voxel") {
          state.activeSubsystems = new Set(["Voxel Engine", "World Gen & Maps"]);
          state.minDegreeFilter = 0;
        } else if (preset === "net") {
          state.activeSubsystems = new Set(["Networking", "Game Loop & Lifecycle"]);
          state.minDegreeFilter = 0;
        } else if (preset === "combat") {
          state.activeSubsystems = new Set(["Player & Combat", "Bot AI"]);
          state.minDegreeFilter = 0;
        }
        document.getElementById("degree-slider").value = state.minDegreeFilter;
        document.getElementById("degree-val").textContent = state.minDegreeFilter;
        applyFilters();
        resetView();
      };
    });

    const degreeSlider = document.getElementById("degree-slider");
    degreeSlider.addEventListener("input", (e) => {
      const val = parseInt(e.target.value, 10);
      state.minDegreeFilter = val;
      document.getElementById("degree-val").textContent = val;
      applyFilters();
    });

    const searchInput = document.getElementById("search-input");
    const autoList = document.getElementById("autocomplete");

    window.addEventListener("keydown", (e) => {
      if (e.key === "/" && document.activeElement !== searchInput) {
        e.preventDefault();
        searchInput.focus();
        searchInput.select();
      }
    });

    searchInput.addEventListener("input", () => {
      const query = searchInput.value.trim().toLowerCase();
      if (!query) {
        autoList.style.display = "none";
        return;
      }
      const matches = RAW_NODES.filter(n => n.name.toLowerCase().includes(query)).slice(0, 10);
      if (matches.length === 0) {
        autoList.style.display = "none";
        return;
      }
      autoList.innerHTML = "";
      matches.forEach(m => {
        const item = document.createElement("div");
        item.className = "autocomplete-item";
        item.innerHTML = `<div><strong style="color:#fff;">${m.name}()</strong><div style="font-size:10px;color:var(--text-muted);">${m.subsystem} • L${m.start_line}</div></div><span class="badge" style="background:${SUBSYSTEMS[m.subsystem]?.badge_bg};color:${SUBSYSTEMS[m.subsystem]?.color};">${m.total_degree} calls</span>`;
        item.onclick = () => {
          searchInput.value = m.name;
          autoList.style.display = "none";
          jumpToNode(m.name);
        };
        autoList.appendChild(item);
      });
      autoList.style.display = "block";
    });

    document.addEventListener("click", (e) => {
      if (!searchInput.contains(e.target) && !autoList.contains(e.target)) {
        autoList.style.display = "none";
      }
    });

    function resetView() {
      state.zoom = 0.85;
      state.panX = width / 2;
      state.panY = height / 2;
    }
    document.getElementById("btn-reset").onclick = resetView;
    document.getElementById("zoom-in").onclick = () => { state.zoom *= 1.25; };
    document.getElementById("zoom-out").onclick = () => { state.zoom *= 0.8; };
    document.getElementById("zoom-fit").onclick = resetView;

    const btnLabels = document.getElementById("btn-labels");
    btnLabels.onclick = () => {
      state.showAllLabels = !state.showAllLabels;
      btnLabels.textContent = state.showAllLabels ? "🏷️ Adaptive Labels" : "🏷️ All Labels";
      btnLabels.classList.toggle("active", state.showAllLabels);
    };

    const btnPhys = document.getElementById("btn-physics-toggle");
    btnPhys.onclick = () => {
      state.physicsRunning = !state.physicsRunning;
      btnPhys.textContent = state.physicsRunning ? "⏸️ Freeze" : "▶️ Unfreeze";
      btnPhys.classList.toggle("active", !state.physicsRunning);
    };

    const flowModal = document.getElementById("flow-modal");
    document.getElementById("btn-flow").onclick = () => {
      renderFlowTable();
      flowModal.style.display = "block";
    };
    document.getElementById("close-flow-modal").onclick = () => {
      flowModal.style.display = "none";
    };

    function renderFlowTable() {
      const table = document.getElementById("flow-table-content");
      const subs = Object.keys(SUBSYSTEMS);
      let html = "<thead><tr><th>Caller \\ Callee</th>";
      subs.forEach(s => {
        html += `<th style="color:${SUBSYSTEMS[s].color};">${s}</th>`;
      });
      html += "</tr></thead><tbody>";

      subs.forEach(s1 => {
        html += `<tr><th style="text-align:left;color:${SUBSYSTEMS[s1].color};">${s1}</th>`;
        subs.forEach(s2 => {
          const cnt = FLOW_MATRIX[s1]?.[s2] || 0;
          const cls = cnt > 0 ? "flow-cell has-calls" : "flow-cell";
          html += `<td class="${cls}">${cnt > 0 ? cnt : "-"}</td>`;
        });
        html += "</tr>";
      });
      html += "</tbody>";
      table.innerHTML = html;
    }

    resizeCanvas();
    initNodesAndPositions();
    setupSubsystemFilters();
    applyFilters();

    const mainNode = state.nodeMap.get("main") || state.nodes[0];
    selectNode(mainNode);

    draw();
  </script>
</body>
</html>
"""

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_template)

    print(f"Generated interactive call graph: {output_path} ({len(html_template)} bytes)")


def main():
    source_file = "fps_ray.c"
    if len(sys.argv) > 1:
        source_file = sys.argv[1]

    if not os.path.isfile(source_file):
        print(f"Error: file not found: {source_file}")
        sys.exit(1)

    print(f"Analyzing {source_file}...")
    functions, cleaned, source = parse_functions(source_file)
    print(f"Discovered {len(functions)} top-level functions.")

    print("Extracting caller-callee dependency edges...")
    edges, flow_matrix = extract_call_graph(functions, cleaned)
    print(f"Resolved {len(edges)} caller-callee call edges.")

    ws_output = "fps_ray_callgraph.html"
    generate_html(functions, edges, flow_matrix, ws_output)

    artifact_dir = "/Users/migr4479/.gemini/antigravity/brain/2c55267e-0b05-4b93-aee3-323c0ad86166"
    if os.path.isdir(artifact_dir):
        artifact_output = os.path.join(artifact_dir, "fps_ray_callgraph.html")
        generate_html(functions, edges, flow_matrix, artifact_output)

    print("Done!")


if __name__ == "__main__":
    main()
