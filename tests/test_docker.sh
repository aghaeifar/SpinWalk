SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd "$SCRIPT_DIR" > /dev/null

docker system prune -a --volumes --force

rm -rf "./build"

echo "Building GPU Docker image..."
docker build  --no-cache -f ../Dockerfile -t spinwalk_gpu ..

rm -rf "./build"
echo "Building CPU Docker image..."
docker build  --no-cache -f ../Dockerfile_CPU -t spinwalk_cpu ..