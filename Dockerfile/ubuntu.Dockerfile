From nvidia/cuda:12.3.1-runtime-ubuntu20.04

RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=America/Los_Angeles apt-get -y install tzdata

RUN apt update && apt install -y wget vim python3-dev python3-numpy g++ build-essential \
	libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libglu1-mesa-dev \
	freeglut3-dev mesa-common-dev libglfw3 libglfw3-dev

RUN wget https://github.com/Kitware/CMake/releases/download/v3.27.9/cmake-3.27.9-linux-x86_64.sh && \
	chmod +x cmake-3.27.9-linux-x86_64.sh && \
	./cmake-3.27.9-linux-x86_64.sh --skip-license --prefix=/usr/local

# https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&Distribution=Ubuntu&target_version=20.04&target_type=deb_local
RUN wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin && \
	mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600 && \
	wget https://developer.download.nvidia.com/compute/cuda/12.3.2/local_installers/cuda-repo-ubuntu2004-12-3-local_12.3.2-545.23.08-1_amd64.deb && \
	dpkg -i cuda-repo-ubuntu2004-12-3-local_12.3.2-545.23.08-1_amd64.deb && \
	cp /var/cuda-repo-ubuntu2004-12-3-local/cuda-*-keyring.gpg /usr/share/keyrings/ && \
	apt-get update && \
	apt-get -y install cuda-toolkit-12-3

WORKDIR /data