Install dependencies:

```
$ sudo apt-get update
$ sudo apt-get install build-essential debhelper git-core pristine-tar git-buildpackage libbz2-dev zlib1g-dev libxml2-dev wget openjdk-6-jdk
```

Install natty dependencies if targetting internal apt repo:

```
$ sudo dpkg -i dpkg_1.16.0~ubuntu7_amd64.deb libalgorithm-diff-perl_1.19.02-1_all.deb libdpkg-perl_1.16.0~ubuntu7_all.deb dpkg-dev_1.16.0~ubuntu7_all.deb libalgorithm-merge-perl_0.08-1_all.deb
```

Create workspace and clone repos. It may be necessary for these repos to exist in adjacent directories.

```
$ export WORKSPACE=/vagrant/workspace
$ mkdir -p $WORKSPACE && cd $WORKSPACE
$ for repo in "ncbi/ncbi-vdb" "ncbi/ngs" "genome-vendor/sra-sdk" ; do git clone https://github.com/$repo ; done
```

Build support libraries. These install to `$HOME/ncbi-outdir` and then are used when building the package. Build order is important.

```
$ for dir in "ngs/ngs-sdk" "ngs/ngs-java" "ncbi-vdb" "ngs" ; do cd $WORKSPACE/$dir && ./configure && make ; done
```

Now you can make sra-sdk (remove the `-uc -us` option if you want to sign during this step)

```
$ cd $WORKSPACE/sra-sdk
$ for repo in "upstream" "pristine-tar" "patch-queue/master" ; do git checkout -b {,origin/}$repo ; done
$ git checkout master
$ git-buildpackage -uc -us --changes-option='-DDistribution=lucid-genome-development'
```

To import a new upstream using gbp:
```
$ git-import-orig --upstream-version=2.4.3 --pristine-tar ../sra-tools-2.4.3.tar.gz
```
