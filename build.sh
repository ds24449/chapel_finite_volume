#!/bin/bash
echo "Compiling finite_volume.chpl"
chpl src/finite_volume.chpl
mkdir -p Generated
mv -f finite_volume Generated/finite_volume
echo "Executable in Generated directory"