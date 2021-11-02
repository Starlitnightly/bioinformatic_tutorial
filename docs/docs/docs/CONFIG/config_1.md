# 配置1：Python环境

对于生物信息学相关的分析，这里分为Python，R语言与linux三部分

## 1. 前言

在python部分，我们主要使用到numpy，pandas，seaborn，matplotlib，scipy，statsmodel这几个包，值得庆幸的是，只需要通过Anaconda的安装，这几个包就可以一并被安装上了，避免各种各样报错的烦恼

## 2. Anaconda下载

如果你没有科学上网的话，建议从清华大学镜像网下载Anaconda3的安装包进行安装，省事儿

下载地址：https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/Anaconda3-5.3.1-Windows-x86_64.exe

## 3. Anaconda安装（转自知乎BG大龍）

### 3.1 安装过程

![img](https://pic1.zhimg.com/80/v2-1a95c6756d90ce6dd74a9f08f6dd50a8_720w.jpg)

![img](https://pic2.zhimg.com/80/v2-29778a46617e491adb554ac5fa5823b1_720w.jpg)

接着就是路径，提醒小白，安装到C盘真的可以避免后续的很多小问题，但是尽管这样我也没有尝试过把它装入C盘。

> 我选择了E盘，单独创建一个文件夹命名为“Anaconda”.
> 注意路径要简单，我的是 E:\Anaconda\ 
>
> ——不要有空格！！！
>
> ——不要有中文字符！！！

![img](https://pic2.zhimg.com/80/v2-64590b21362f65132f54d1a597f0c809_720w.jpg)

大家注意到这一幅图跟上一幅图不一样，那是因为，我建议两个都打勾，第一个不打勾对初学者其实挺不友好的，打勾了能省很多事儿

![img](https://pic2.zhimg.com/80/v2-8c295d70c79db8487180a10a3b1798a1_720w.jpg)

点击install，等待不太漫长的进度条……

![img](https://pic1.zhimg.com/80/v2-76ba09f793eb5ed5044e63a9a07bb4c0_720w.jpg)

提示安装成功……

![img](https://pic4.zhimg.com/80/v2-fde78f1bb363c59f31ed653e872e788b_720w.jpg)

提示安装VScode，选择点击“skip”

![img](https://pic1.zhimg.com/80/v2-ec381c321d4a6abc833acc39e2ab1698_720w.jpg)

两个“learn”，都取消打勾

![img](https://pic2.zhimg.com/80/v2-6ed983e8cedf48dc0c0870d3de2c620d_720w.jpg)

### 3.2 检验安装

我们在cmd中输入python，检查是否有Python环境

> Q：cmd怎么打开
>
> A：“cmd命令的打开方法：1、在电脑桌面中使用“WIN＋R”组合键，打开的“运行”窗口，输入“cmd”命令并回车即可打开；2、打开“开始”菜单，在搜索框中输入“cmd”，点击“cmd.exe”即可打开。”

![img](https://pic3.zhimg.com/80/v2-bad7f81e99c68cce4359359f40bcaa0a_720w.jpg)

在cmd中输入： conda info，——查看是否有？ (检验安装成功的标志)

![img](https://pic3.zhimg.com/80/v2-26f33b7c61cc29739f172bb1a1392d66_720w.jpg)

### 3.3 检验其他是否安装成功。尤其是 Anaconda Navifator

> 点击，看是否能够进入界面，若成功，大功告成。

![img](https://pic1.zhimg.com/80/v2-724a6e563be81c17c9cd060ddc8e8d50_720w.jpg)

![img](https://pic1.zhimg.com/80/v2-4b2b25ef579b8dfdd55e5f4950de18cc_720w.jpg)



到这里，你已经完成了Python环境的基础配置，至于更高阶的配置，那么我们等到需要进行对应的分析的时候，再进行安装。