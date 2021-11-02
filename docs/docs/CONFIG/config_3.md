# 配置3：R语言环境

一方面，由于我的个人习惯，我喜欢用jupyter来处理代码块；另一方面，作者给的代码全是R语言，并且都是在服务器上运行，所以我们需要一个linux环境，这里，我用了windows的wsl，所以本文档主要是讲，wsl如何安装jupyter，以及如何添加IRkernel，并且装上所有的依赖（这个真的很难，各种报错）

## 1、wsl下载

（1）打开 Microsoft Store , 搜索 Linux 会有三个结果 Ubuntu , openSUSE Leap42 , SUSE Linux Enterprise Srever

![img](https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/clip_image002.png)

（2）选择ubuntu后一键安装直到完成，在通知里打开即可

 

![img](https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/clip_image004.jpg)

（3）到 控制面板\所有控制面板项\程序和功能 中选择 启用或者关闭Windows功能

![img](https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/clip_image006.png)

（4）找到 适用于Linux的Windows子系统 并勾选,然后 确定 选择 立即重新启动

![img](https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/clip_image007.png)

（5）重启之后再次打开 ***Ubuntu*** 就可以使用了,根据提示输入用户名,两次输入密码就可以看见熟悉的命令行

![img](https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/clip_image009.jpg)

## 2、下载源设置

由于linux下大部分安装包会直接指向国外的下载地址，所以我们一般会更改下载源，改成国内的，你可以理解成，国内一些大厂以及高校都将linux的包全部下好了，我们去他那里下载

### 2.1 编辑

首先我们要打开下载源编辑器

```
 vim /etc/apt/sources.list
```

打开后需要清空该下载源

```
#输入%d后直接按回车
%d
```

然后黏贴新的下载源

```python
deb https://mirrors.aliyun.com/ubuntu/ bionic main restricted universe multiverse
deb https://mirrors.aliyun.com/ubuntu/ bionic-security main restricted universe multiverse
deb https://mirrors.aliyun.com/ubuntu/ bionic-updates main restricted universe multiverse
deb https://mirrors.aliyun.com/ubuntu/ bionic-proposed main restricted universe multiverse
deb https://mirrors.aliyun.com/ubuntu/ bionic-backports main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/linux/ubuntu/ bionic-cran40/
```

（4）保存退出

```
#按下esc后输入:wq,然后回车
:wq
```

### 2.2 更新下载

```
sudo apt-get update
```

### 2.3 公钥设置（可选）

理论上，（5）的命令你输入了之后可能会报错，报错之后我们输入下面这行代码就好了

```
apt install ca-certificates
curl -sL "http://keyserver.ubuntu.com/pks/lookup?search=0x51716619E084DAB9&op=get" | sudo apt-key add
#上面这行代码输入完后，再输入一次下面这个
sudo apt-get update
```



## 2、代理设置（可选）

怎么说呢，我是wsl，然后电脑访问github得科学上网，于是我的代理端口是4780，但是wsl不会自动识别啊，所以我们得把端口映射过去

### 2.1 安装

我们这里安装的是proxychains4

```
sudo apt update
sudo apt install proxychains4
```

### 2.2 配置

```
sudo vim /etc/proxychains.conf
```

注释掉socks4 127.0.0.1那一行，在最后加上代理工具的设置，如：

```
socks4  127.0.0.1 4780
```

然后ping www.github.com 就能成功了

## 3、miniconda安装

虽然说吧，wsl不支持conda的虚拟环境，但是我主要是想要一个python环境，于是就安装了miniconda。

### 3.1 新建目录

首先我们新建一个文件夹，命令行输入这两行命令

```python
mkdir anaconda
cd anaconda
```

### 3.2 下载miniconda

然后从清华的镜像下载miniconda，命令行输入wget

```
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh
```

### 3.3 安装

下载好后，我们直接运行安装程序，有点像windows下的双击安装，只不过这里是以命令的形式

```
bash Miniconda3-4.7.12.1-Linux-x86_64.sh
```

等一会儿就安装好了

## 4、jupyter lab安装

我这里讲的是wsl下怎么安装，所以就跟服务器上的有一些区别，但其实也不大，服务器上我更喜欢用conda来装东西

### 4.1 安装

直接使用pip命令安装，然后等待安装完成

```
sudo pip install jupyterlab
```

### 4.2 生成配置

我们首先生成一个notebook的配置文件

```
jupyter notebook --generate-config --allow-root
```

上述代码执行成功后，会出现如下信息

```
#Writing default config to: /root/.jupyter/jupyter_notebook_config.py
```

### 4.3 生成密码

在jupyter 5.0以后，我们通过浏览器访问笔记本都需要输入密码，所以我们这里设置一下

打开ipython

```
ipython
```

执行下面的内容，***就是你的密码，我这里输入的是123

```python
In [1]: from notebook.auth import passwd
In [2]: passwd()
Enter password:***
Verify password:***
Out[2]: 'sha1:67c9e60bb8b6:9ffede0825894254b2e042ea597d771089e11aed'
```

`sha1:67c9e60bb8b6:9ffede0825894254b2e042ea597d771089e11aed` 这一串就是要在 `jupyter_notebook_config.py` 添加的密码。

### 4.4 更新配置

我们在上面的步骤中已经生成了配置文件与密码，接下来就是将其导入配置的过程，在命令行输入下面的命令，打开配置

```
vim /root/.jupyter/jupyter_notebook_config.py
```

在 `jupyter_notebook_config.py` 中找到下面的行，取消注释并修改。

```
c.NotebookApp.ip='*'#163行
c.NotebookApp.password = u'sha:ce...刚才复制的那个密文'  #217行
c.NotebookApp.open_browser = False#208
c.NotebookApp.port =8888 #可自行指定一个端口, 访问时使用该端口228行
```

以及密码行

```
c.NotebookApp.password = u'sha1:67c9e60bb8b6:9ffede0825894254b2e042ea597d771089e11aed'
```

设置完之后就是常规的vim保存退出

### 4.5 验证jupyter-lab的安装

我们随便前进一个目录，注意不能是在隐藏目录 (以 . 开头的目录)下启动 jupyter notebook, 否则无法正常访问文件。输入

```
jupyter lab
```

命令提示会出现下面的命令

```
[I 2021-10-06 22:41:13.633 LabApp] JupyterLab application directory is /root/miniconda3/share/jupyter/lab
[I 2021-10-06 22:41:13.637 ServerApp] jupyterlab | extension was successfully loaded.
[I 2021-10-06 22:41:13.638 ServerApp] The port 8888 is already in use, trying another port.
[I 2021-10-06 22:41:13.639 ServerApp] Serving notebooks from local directory: /mnt/c/Windows/System32
[I 2021-10-06 22:41:13.639 ServerApp] Jupyter Server 1.11.1 is running at:
[I 2021-10-06 22:41:13.639 ServerApp] http://localhost:8889/lab
[I 2021-10-06 22:41:13.640 ServerApp]  or http://127.0.0.1:8889/lab
[I 2021-10-06 22:41:13.640 ServerApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation)
```

我们在浏览器输入http://localhost:8889/lab 就能打开我们的网址，输入之前设置的密码123，就可以正常使用了

## 5、IRkernel安装

我们目的是在jupyter里面使用R语言，所以需要在jupyter里面配置R环境

### 5.1 R安装

我在前面的下载源里面就设置好了R4.0以上版本的安装CRAN40，所以我们直接apt-get就能下载到最新的R包了

```
sudo apt-get install r-base
```

### 5.2 安装lib

由于这个包需要一定的lib依赖，所以我们简单配置一下

```
sudo apt-get install libzmq3-dev libssl-dev  libcurl4-openssl-dev libxml2-dev
#安装完后进入R环境
sudo R
```

在R环境里面，我们先输入如下语句

```
Sys.setenv(R_INSTALL_STAGED = FALSE)
```

就可以避免一些bug，重命名bug

```
install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'devtools', 'uuid', 'digest', 'psych'))
```

等上述的包安装完之后，输入

```R
devtools::install_github('IRkernel/IRkernel')

IRkernel::installspec()
```

如果显示如下结果，说明配置成功。 

`[InstallKernelSpec] Installed kernelspec ir in /home/xxx/.local/share/jupyter/kernels/ir`