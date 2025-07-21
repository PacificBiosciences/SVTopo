# Docker build
This build script uses a modified version of the official Rust Docker image to compile this crate.
If enabled, it will build a static binary that can be uploaded or shared as necessary.
This system is meant to be relatively generic such that it can be easily copy/pasted from one crate to the next.
This requires configuring the `TOOL_NAME` in Bamboo or if you are running a local install.

There is both a Bamboo and local build system in this folder:
- `Dockerfile` - used by both Bamboo and local to create the Docker image
- `bamboo_docker_build.sh` - run on Bamboo inside the Docker container, this is the workhorse for testing and building
- `build_image.sh` - shorcut script for building the image locally
- `local_x86_64-unknown-linux-gnu.sh` - shortcut script for creating a release locally, this does NOT do testing

Steps to build locally:
```bash
# ssh to your resource with docker installed
# configure tool name
export TOOL_NAME="{name}"

# one-time setup to create build image if new Docker instance
./build_image.sh

# run the build script
cd ${ROOT}/docker-rust
./local_x86_64-unknown-linux-gnu.sh

# check output help works
./release/${TOOL_NAME}-v${VERSION}-x86_64-unknown-linux-gnu/${TOOL_NAME} -h
```
