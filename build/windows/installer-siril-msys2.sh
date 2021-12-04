# Install Inno Setup.
#wget https://jrsoftware.org/download.php/is.exe
#./is.exe //SILENT //SUPPRESSMSGBOXES //CURRENTUSER //SP- //LOG="innosetup.log"

# Construct now the installer.
MAJOR_VERSION=`grep "version : " meson.build | head -n 1  | awk '{ print $3 }' | tr -d "'," | cut -d . -f 1`
MINOR_VERSION=`grep "version : " meson.build | head -n 1  | awk '{ print $3 }' | tr -d "'," | cut -d . -f 2`
MICRO_VERSION=`grep "version : " meson.build | head -n 1  | awk '{ print $3 }' | tr -d "'," | cut -d . -f 3`
cd build/windows/installer
./compile.bat ${MAJOR_VERSION}.${MINOR_VERSION}.${MICRO_VERSION} ../../.. siril-w64 ../../.. siril-w64

# Test if the installer was created and return success/failure.
if [ -f "_Output/siril-${MAJOR_VERSION}.${MINOR_VERSION}.${MICRO_VERSION}-setup.exe" ]; then
  cd _Output/
  INSTALLER="siril-${MAJOR_VERSION}.${MINOR_VERSION}.${MICRO_VERSION}-setup.exe"
  sha256sum $INSTALLER > ${INSTALLER}.SHA256SUMS
  sha512sum $INSTALLER > ${INSTALLER}.SHA512SUMS
  exit 0
else
  exit 1
fi
