编译镜像的dockerfile在Dockerfile文件夹下，按照windows wsl2和ubuntu环境分为wsl.Dockerfile和ubuntu.Dockerflie。 使用方法如下：

1. 安装nvidia-docker，可参考如下文档：https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html
```
apt-get update
apt-get install -y nvidia-container-toolkit
```

2. 编译镜像(以windows wsl2环境为例)
```
docker build -f Dockerfile/wsl.Dockerfile -t {your-image-name(lower case)} .
```

3. 启动镜像
```
docker run -it --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix --gpus all --net=host -v {Path of GPU-PBD and TinyVisualizer(for example /PBD)}:/data:rw {your-image-name} /bin/bash
```
启动后，即可在该container中编译以及运行GPU-PBD和TinyVisualizer。
