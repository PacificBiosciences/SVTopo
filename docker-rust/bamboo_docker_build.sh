#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail

# binary CLI should have "true"; crate/libraries should have "false"
BINARY_RELEASE=true

# make sure the tool name is set externally
if ! [ -n "${TOOL_NAME}" ]; then
    echo "TOOL_NAME is not set"
    exit 1;
fi 

# these components do not need to run in the static configuration
echo "cargo clippy"
cargo clippy -- -D warnings

# doc tests do not behave inside static mode, nor do they need to since it's for documentation
echo "cargo test --doc"
cargo test --release --doc

# license checker test
echo "checking license"

# check that our license file is up to date
cargo license -j | \
    /license_checker/scripts/check_licenses.py
echo "license checks passed."

# run the functional tests
echo "functional tests"
eval "$(/usr/local/bin/micromamba shell hook --shell bash)"
micromamba activate py310
bash test/scripts/run_end_to_end_tests.sh

# everything after here we want to run static targeting x86_64-unknown-linux-gnu
export RUSTFLAGS="-C target-feature=+crt-static -C relocation-model=static"
echo "cargo test"
cargo test --release --lib --target x86_64-unknown-linux-gnu -- --include-ignored

if $BINARY_RELEASE; then
    # we want to build the full release binary
    echo "cargo build"
    cargo build --release --target x86_64-unknown-linux-gnu

    # make and populate the releases
    RELEASE_FOLDER="release"
    mkdir -p ${RELEASE_FOLDER}

    # make and populate this particular version
    BINARY_FILE=./target/x86_64-unknown-linux-gnu/release/${TOOL_NAME}
    VERSION=`${BINARY_FILE} -V | awk -F '[ -]' '{print $2}'`
    TAR_FOLDER="${TOOL_NAME}-v${VERSION}-x86_64-unknown-linux-gnu"
    TARBALL="${TAR_FOLDER}.tar.gz"
    VERSION_FOLDER="${RELEASE_FOLDER}/${TAR_FOLDER}"
    mkdir -p ${VERSION_FOLDER}
    cp ${BINARY_FILE} ${VERSION_FOLDER}/${TOOL_NAME}

    # generate an md5sum
    (cd ${VERSION_FOLDER} && \
        md5sum ${TOOL_NAME} > ${TOOL_NAME}.md5)

    # finally make the tarball
    (cd ${RELEASE_FOLDER} && \
        tar -czvf ${TARBALL} ${TAR_FOLDER})
fi

