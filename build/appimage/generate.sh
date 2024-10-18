#!/bin/bash
set -e

# Create application root directory
APP_DIR="MyApp.AppDir"
mkdir -p "$APP_DIR"
cd "$APP_DIR"

# Create directory structure
mkdir -p usr/bin
mkdir -p usr/lib
mkdir -p usr/lib/python3.10
mkdir -p usr/lib/python3.10/site-packages

# Copy Python executable and libraries
PYTHON_VERSION="3.10"  # Adjust this to match your Python version
cp /usr/bin/python${PYTHON_VERSION} usr/bin/
cp /usr/bin/python3 usr/bin/
cp /usr/bin/python usr/bin/

# Copy Python standard library
cp -r /usr/lib/python${PYTHON_VERSION}/* usr/lib/python${PYTHON_VERSION}/

# Copy required shared libraries
ldd usr/bin/python${PYTHON_VERSION} | grep "=> /" | awk '{print $3}' | xargs -I '{}' cp -v '{}' usr/lib/

# Copy additional Python dependencies
cp -r /usr/lib/python3/dist-packages usr/lib/python${PYTHON_VERSION}/
cp -r /usr/lib/python${PYTHON_VERSION}/site-packages/* usr/lib/python${PYTHON_VERSION}/site-packages/ || true

# Copy your AppRun script
cat > AppRun << 'EOL'
#!/bin/bash
HERE="$(dirname "$(readlink -f "${0}")")"

# Set up Python environment
export PYTHONHOME="${HERE}/usr"
export PYTHONPATH="${HERE}/usr/lib/python3.10:${HERE}/usr/lib/python3.10/site-packages:${PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1
export LD_LIBRARY_PATH="${HERE}/usr/lib:${LD_LIBRARY_PATH}"

# Set up other environment variables
export PATH="${HERE}/usr/bin:${PATH}"
export XDG_DATA_DIRS="${HERE}/usr/share:${XDG_DATA_DIRS:-/usr/local/share:/usr/share}"

# Execute the main application
if [ ! -z "$1" ] && [ -e "$HERE/usr/bin/$1" ] ; then
    MAIN="$HERE/usr/bin/$1"
    shift
else
    MAIN="$HERE/usr/bin/python3"  # Default to python3 if no specific command
fi

exec "$MAIN" "$@"
EOL

chmod +x AppRun

# Create .desktop file
cat > MyApp.desktop << 'EOL'
[Desktop Entry]
Name=MyApp
Exec=python3
Icon=python
Type=Application
Categories=Development;
EOL

# Copy Python icon
cp /usr/share/pixmaps/python3.xpm python.xpm

# Fix permissions
chmod +x AppRun
chmod +x usr/bin/*

echo "AppDir structure created. Now you can create the AppImage using appimagetool:"
echo "appimagetool MyApp.AppDir MyApp-x86_64.AppImage"
