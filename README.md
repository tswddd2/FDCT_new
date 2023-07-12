# FDCT_new

## Environment Setup (Linux)

Download and unzip boost 1.55 from this [link](https://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.gz/download)

```
cd ${boost_path}
wget -O boost_1_55_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.gz/download
tar xzvf boost_1_55_0.tar.gz
```

## run experiment
```
cd ${FDCT_path}
python3 exp.py ${boost_path}/boost_1_55_0
```
