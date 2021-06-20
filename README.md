# 2021 春季算法期末项目

任务描述：[PROJECT.pdf](./PROJECT.pdf)

算法描述：[ARCHITECTURE](./ARCHITECTURE.md)

## 依赖项

* `*`：
    * CMake 3.20.3
    * GNU Make 4.3
    * GNU G++ 11.1.0 (`-std=c++20`)
    * `-lpthread`
* `locate-demo`：
    * OpenGL 3.3
    * SDL 2.0.14
    * GLEW 2.2.0

## 编译

```shell
mkdir build
cd build
cmake ..
make -j
```

## 运行

### 样例

```shell
cd build
./locate -r ../data/sample/ref.fasta -l ../data/sample/long.fasta -j8 2> sample.locate.txt
./dump -r ../data/sample/ref.fasta -l ../data/sample/long.fasta -p sample.locate.txt -m 1.0 -j8 2> sample.dump.txt
./analyze -r ../data/sample/ref.fasta -l ../data/sample/long.fasta -p sample.locate.txt -d sample.dump.txt 2> sample.answer.txt
```

### 正式数据

```shell
cd build
./locate -r ../data/final/ref.fasta -l ../data/final/long.fasta -j8 2> final.locate.txt
./dump -r ../data/final/ref.fasta -l ../data/final/long.fasta -p final.locate.txt -m 1.0 -j8 2> final.dump.txt
./analyze -r ../data/final/ref.fasta -l ../data/final/long.fasta -p final.locate.txt -d final.dump.txt 2> final.answer.txt
cd ..
ln -s build/final.answer.txt sv.bed
```
