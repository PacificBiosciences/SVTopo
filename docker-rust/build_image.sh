
# make sure the tool name is set externally
if ! [ -n "${TOOL_NAME}" ]; then
    echo "TOOL_NAME is not set"
    exit 1;
fi 

# very simple build script, saves it to "${TOOL_NAME}_build:latest"
version="1.0.0"
docker build . \
    --tag ${TOOL_NAME}_build:${version} \
    --tag ${TOOL_NAME}_build:latest 