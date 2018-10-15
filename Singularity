Bootstrap:docker
From:centos:centos7.4.1708

%help
    Container for nanopore tools

%labels
    MAINTAINER Bruce Moran

%environment
    #environment variables defined in post section using
    ##SINGULARITY_ENVIRONMENT variable

%post
    #essential utilities
    yum -y install git wget bzip2 unzip which emacs

    #language and libraries
    yum -y install java-1.8.0-openjdk-devel gcc gcc-c++ glibc-devel make ncurses ncurses-devel zlib-devel libbzip2-devel bzip2-devel xz-devel perl-DBI perl-core lapack-devel atlas-devel freetype freetype-devel libpng-devel readline-devel pcre-devel libtool openssl-devel libxml2-devel mysql-devel tcl-devel tk-devel readline readline-devel pcre pcre-devel libcurl libcurl-devel

    #libclas and libatlas aren't put in the right places
    ln -sf /usr/lib64/atlas/libtatlas.so /usr/lib64/libatlas.so
    ln -sf /usr/lib64/atlas/libsatlas.so /usr/lib64/libcblas.so

    ##python 3.6
    yum -y install https://centos7.iuscommunity.org/ius-release.rpm
    yum -y install python36u python36u-pip python36u-devel
    pip3.6 install --upgrade pip
    pip3.6 install --upgrade setuptools

    mkdir -p /usr/local/src
    cd /usr/local/src

    ##define env vars via S..._E... env var when in post
    ##see: https://www.sylabs.io/guides/2.5/user-guide/environment_and_metadata.html
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib' >> $SINGULARITY_ENVIRONMENT

    ##setting more that LANG locale is an issue for several tools
    ##https://github.com/CentOS/sig-cloud-instance-images/issues/71
    localedef -i en_US -f UTF-8 en_US.UTF-8
    echo -e "LANGUAGE="C"\nLC_ALL="en_US.utf8"" >> /etc/locale.conf
    echo 'export LANG=en_US.UTF-8' >> $SINGULARITY_ENVIRONMENT
    echo 'export LANGUAGE=C' >> $SINGULARITY_ENVIRONMENT
    echo 'export LC_ALL=C' >> $SINGULARITY_ENVIRONMENT
    echo 'export LC_CTYPE=C' >> $SINGULARITY_ENVIRONMENT

    ##R
    #source
    wget https://cran.rstudio.com/src/base/R-3/R-3.5.1.tar.gz
    tar xf R-3.5.1.tar.gz
    cd R-3.5.1
    ./configure --with-x=no --prefix=/usr/local/
    make
    make install

    #packages
    R --slave -e 'install.packages("BiocManager", repos="https://cloud.r-project.org/")'
    R --slave -e 'library("BiocManager"); lapply(c("QDNAseq", "QDNAseq.hg19", "PSCBS", "matrixStats", "ggplot2", "caTools", "dplyr", "data.table", "futile.logger", "optparse", "plyr", "readr", "reshape2", "scales", "viridis", "yaml"), function(f){BiocManager::install(f, dependencies=TRUE, update=FALSE, ask=FALSE)})'

    cd /usr/local/src

    #albacore
    pip3.6 install https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.3.3-cp36-cp36m-manylinux1_x86_64.whl

    #minimap2
    git clone https://github.com/lh3/minimap2
    cd minimap2 && make
    ln -s /usr/local/src/minimap2/minimap2 /usr/local/bin/minimap2
    cd /usr/local/src

    #MinionQC
    wget https://github.com/roblanf/minion_qc/blob/feature/detect_png/MinIONQC.R -O MinIONQC.R
    echo -e "#! /bin/bash\nMINIONQCR=\$(echo \$0".R")\nRscript --vanilla \$MINIONQCR -i \$1\n" > /usr/local/bin/MinIONQC
    chmod 755 /usr/local/bin/MinIONQC

    #NanoStat
    pip3.6 install nanostat

    #samtools
    wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
    tar xf samtools-1.8.tar.bz2
    cd samtools-1.8
    ./configure --prefix=/usr/local/
    make
    make install
    cd /usr/local/src

    #bcftools
    wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
    tar xf bcftools-1.8.tar.bz2
    cd bcftools-1.8
    ./configure --prefix=/usr/local
    make
    make install
    cd /usr/local/src

    #htslib
    wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
    tar xf htslib-1.8.tar.bz2
    cd htslib-1.8
    ./configure --prefix=/usr/local
    make
    make install
    cd /usr/local/src

    #NextFlow
    curl -s https://get.nextflow.io | bash
    mv nextflow /usr/local/bin/
    chmod 777 /usr/local/bin/nextflow

%runscript
    #if need to run stuff, put here
