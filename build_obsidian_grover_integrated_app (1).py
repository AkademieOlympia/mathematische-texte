#!/usr/bin/env python3
"""
macOS .app Bundle Creator for Obsidian Grover Integrated Studio
----------------------------------------------------------------
Erstellt ein eigenständiges anklickbares macOS Application Bundle für den Mac.
"""

import os
import sys
import stat
import shutil
from pathlib import Path

def build_macos_app():
    app_name = "ObsidianGroverIntegratedStudio"
    build_dir = Path("/workspace/scratch")
    app_dir = build_dir / f"{app_name}.app"
    
    # Bundle Directory Structure
    contents_dir = app_dir / "Contents"
    macos_dir = contents_dir / "MacOS"
    resources_dir = contents_dir / "Resources"
    
    if app_dir.exists():
        shutil.rmtree(app_dir)
        
    os.makedirs(macos_dir, exist_ok=True)
    os.makedirs(resources_dir, exist_ok=True)
    
    # 1. Copy the main streamlit script into Resources
    script_source = build_dir / "obsidian_grover_integrated_studio.py"
    script_target = resources_dir / "obsidian_grover_integrated_studio.py"
    if script_source.exists():
        shutil.copy(script_source, script_target)
    else:
        # Fallback if building directly from artifacts
        art_source = Path("/workspace/artifacts/obsidian_grover_integrated_studio.py")
        if art_source.exists():
            shutil.copy(art_source, script_target)

    # 2. Info.plist
    info_plist = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleExecutable</key>
    <string>launch_app</string>
    <key>CFBundleIdentifier</key>
    <string>com.quantum.obsidian.grover.integrated</string>
    <key>CFBundleName</key>
    <string>{app_name}</string>
    <key>CFBundlePackageType</key>
    <string>APPL</string>
    <key>CFBundleShortVersionString</key>
    <string>1.0</string>
    <key>LSMinimumSystemVersion</key>
    <string>12.0</string>
    <key>NSHighResolutionCapable</key>
    <true/>
</dict>
</plist>
"""
    with open(contents_dir / "Info.plist", "w", encoding="utf-8") as f:
        f.write(info_plist)

    # 3. Launch Script Executable
    launch_script = f"""#!/bin/bash
DIR="$( cd "$( dirname "${{BASH_SOURCE[0]}}" )" && pwd )"
RESOURCES="$DIR/../Resources"

# Native Python / Streamlit Launcher
if command -v streamlit &> /dev/null
then
    streamlit run "$RESOURCES/obsidian_grover_integrated_studio.py" --server.port 8501
else
    python3 -m streamlit run "$RESOURCES/obsidian_grover_integrated_studio.py" --server.port 8501
fi
"""
    executable_path = macos_dir / "launch_app"
    with open(executable_path, "w", encoding="utf-8") as f:
        f.write(launch_script)
        
    # Make launcher executable
    st_mode = os.stat(executable_path).st_mode
    os.chmod(executable_path, st_mode | stat.S_IEXEC)
    
    print(f"macOS .app Bundle erfolgreich erstellt: {app_dir}")
    return app_dir

if __name__ == "__main__":
    build_macos_app()
