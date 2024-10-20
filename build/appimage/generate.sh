# Section 7: AppImage Finalization
########################################################################
# Create the artifacts directory where CI expects to find the files
mkdir -p "$BASE_DIR/build/appimage"

# Find the generated AppImage files
APPIMAGE_FILES=(Siril*.AppImage*)
if [ ${#APPIMAGE_FILES[@]} -eq 0 ]; then
    echo "ERROR: No AppImage files found"
    ls -la
    exit 1
fi

# Move each AppImage file individually
for file in "${APPIMAGE_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "Moving $file to $BASE_DIR/build/appimage/"
        cp "$file" "$BASE_DIR/build/appimage/" || {
            echo "Failed to copy $file to artifacts directory"
            exit 1
        }
        # Verify the copy succeeded before removing the original
        if [ -f "$BASE_DIR/build/appimage/$(basename "$file")" ]; then
            rm "$file"
        else
            echo "ERROR: Failed to verify copied file: $BASE_DIR/build/appimage/$(basename "$file")"
            exit 1
        fi
    else
        echo "WARNING: AppImage file not found: $file"
    fi
done

# Verify the files were moved successfully
if [ "$(ls -A "$BASE_DIR/build/appimage/")" ]; then
    echo "AppImage files successfully moved to artifacts directory:"
    ls -la "$BASE_DIR/build/appimage/"
    echo "AppImage creation and deployment completed successfully"
else
    echo "ERROR: No files found in artifacts directory after move operation"
    exit 1
fi
