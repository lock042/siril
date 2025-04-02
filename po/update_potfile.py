#!/usr/bin/env python3
"""
Cross-platform script to update the POTFILES.in file for Siril.
"""

import os
import re
from pathlib import Path
import time

def search_translation_strings(src_dir, exclude_dirs):
    """Native Python implementation for searching translation strings _('...')"""
    matched_files = []
    total_files = 0
    
    print("Searching for files containing translation strings _('...')")
    start_time = time.time()
    
    # List of extensions to search
    extensions = ['.c', '.cpp']
    
    for ext in extensions:
        for filepath in Path(src_dir).glob(f'**/*{ext}'):
            total_files += 1
            path_str = str(filepath)
            
            # Check if path contains an excluded directory
            if any(exclude_dir in path_str for exclude_dir in exclude_dirs):
                continue
                
            # Check file content
            try:
                with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                    # Search for translation patterns (_("..."))
                    content = f.read()
                    # Regex that searches for _(" followed by ") with anything in between
                    # This regex supports multi-line and concatenated strings
                    pattern = r'_\s*\(\s*".*?"\s*\)'
                    if re.search(pattern, content, re.DOTALL):
                        rel_path = os.path.relpath(filepath).replace('\\', '/')
                        matched_files.append(rel_path)
            except Exception as e:
                print(f"Error reading {filepath}: {e}")
    
    elapsed_time = time.time() - start_time
    print(f"Search completed in {elapsed_time:.2f} seconds, {total_files} files analyzed")
    
    return matched_files

def main():
    # Configuration
    potfiles_path = 'po/POTFILES.in'
    src_dir = 'src'
    exclude_dirs = ['tests', 'build', 'lcms_acceleration', 'deconvolution', 'avi_pipp', 'registration/matching']
    
    print(f"Updating {potfiles_path}...")
    
    # Read current content of POTFILES.in
    with open(potfiles_path, 'r', encoding='utf-8') as f:
        current_content = f.read()
    
    # Extract header (comments and instructions)
    header_match = re.search(r'^(.*?# after adding files here\.).*', current_content, re.DOTALL | re.MULTILINE)
    if not header_match:
        print("Warning: Unrecognized header format in POTFILES.in")
        header = "# List of source files containing translatable strings.\n# Please keep this file in alphabetical order; run ./sort-potfiles\n# after adding files here."
    else:
        header = header_match.group(1)
    
    # Extract UI section
    ui_section_match = re.search(r'# the next lines of ui MUST always be here.*?\n[^\n]*\.ui(\n[^\n]*\.ui)*', current_content, re.DOTALL)
    if not ui_section_match:
        print("Warning: UI section not found")
        ui_section = "# the next lines of ui MUST always be here (it is not sorted)"
    else:
        ui_section = ui_section_match.group(0)
    
    print(f"UI section extracted: {ui_section.count('.ui')} UI files")
    
    try:
        # Use native Python method to search for translation strings
        matched_files = search_translation_strings(src_dir, exclude_dirs)
        
        # Sort files
        matched_files.sort()
        
        print(f"Found {len(matched_files)} source files with translation strings")
        
        # Display first few lines for verification
        if matched_files:
            print("\nFirst files found (5 first):")
            for f in matched_files[:5]:
                print(f"  - {f}")
            
            if len(matched_files) > 5:
                print(f"  ... and {len(matched_files) - 5} other files")
        
        # Create new file content
        new_content = f"{header}\n{ui_section}\n" + '\n'.join(matched_files)
        
        # Write new content to file
        with open(potfiles_path, 'w', encoding='utf-8') as f:
            f.write(new_content)
        
        print(f"\nFile {potfiles_path} updated successfully.")
        print(f"  - UI section: {ui_section.count('.ui')} files")
        print(f"  - Source section: {len(matched_files)} files")
    
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
