#!/bin/sh

if [ $# -ne 2 ]
then
    echo "Usage: copydirsonly.sh <source_dir> <target_dir>"
    exit 1
fi

rsync -av -f"+ */" -f"- *" "${1}" "${2}"
exit 0
