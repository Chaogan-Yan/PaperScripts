# Use Ubuntu 16.04 LTS
FROM ubuntu:xenial-20200114

MAINTAINER Chao-Gan Yan <ycg.yan@gmail.com>

# Pre-cache neurodebian key
# COPY docker/files/neurodebian.gpg /usr/local/etc/neurodebian.gpg

# Prepare environment
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    curl \
                    bzip2 \
                    ca-certificates \
                    xvfb \
                    build-essential \
                    autoconf \
                    libtool \
                    pkg-config \
                    git && \
    curl -sL https://deb.nodesource.com/setup_10.x | bash - && \
    apt-get install -y --no-install-recommends \
                    nodejs && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# Installing freesurfer
RUN curl -sSL https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.1/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.1.tar.gz | tar zxv --no-same-owner -C /opt \
    --exclude='freesurfer/diffusion' \
    --exclude='freesurfer/docs' \
    --exclude='freesurfer/fsfast' \
    --exclude='freesurfer/lib/cuda' \
    --exclude='freesurfer/lib/qt' \
    --exclude='freesurfer/matlab' \
    --exclude='freesurfer/mni/share/man' \
    --exclude='freesurfer/subjects/fsaverage_sym' \
    --exclude='freesurfer/subjects/fsaverage3' \
    --exclude='freesurfer/subjects/fsaverage4' \
    --exclude='freesurfer/subjects/cvs_avg35' \
    --exclude='freesurfer/subjects/cvs_avg35_inMNI152' \
    --exclude='freesurfer/subjects/bert' \
    --exclude='freesurfer/subjects/lh.EC_average' \
    --exclude='freesurfer/subjects/rh.EC_average' \
    --exclude='freesurfer/subjects/sample-*.mgz' \
    --exclude='freesurfer/subjects/V1_average' \
    --exclude='freesurfer/trctrain'

ENV FSL_DIR="/usr/share/fsl/5.0" \
    OS="Linux" \
    FS_OVERRIDE=0 \
    FIX_VERTEX_AREA="" \
    FSF_OUTPUT_FORMAT="nii.gz" \
    FREESURFER_HOME="/opt/freesurfer"
ENV SUBJECTS_DIR="$FREESURFER_HOME/subjects" \
    FUNCTIONALS_DIR="$FREESURFER_HOME/sessions" \
    MNI_DIR="$FREESURFER_HOME/mni" \
    LOCAL_DIR="$FREESURFER_HOME/local" \
    MINC_BIN_DIR="$FREESURFER_HOME/mni/bin" \
    MINC_LIB_DIR="$FREESURFER_HOME/mni/lib" \
    MNI_DATAPATH="$FREESURFER_HOME/mni/data"
ENV PERL5LIB="$MINC_LIB_DIR/perl5/5.8.5" \
    MNI_PERL5LIB="$MINC_LIB_DIR/perl5/5.8.5" \
    PATH="$FREESURFER_HOME/bin:$FSFAST_HOME/bin:$FREESURFER_HOME/tktools:$MINC_BIN_DIR:$PATH"


# Set CPATH for packages relying on compiled libs (e.g. indexed_gzip)
ENV PATH="/usr/local/miniconda/bin:$PATH" \
    CPATH="/usr/local/miniconda/include/:$CPATH" \
    LANG="C.UTF-8" \
    LC_ALL="C.UTF-8" \
    PYTHONNOUSERSITE=1



# Update system and istall pakages
RUN apt-get -qq update && apt-get -qq install -y \
    x11vnc xvfb suckless-tools stterm connectome-workbench parallel wget unzip && \
    apt-get update


    
# Install MATLAB MCR
ENV MATLAB_VERSION R2020a
RUN mkdir /opt/mcr_install && \
    mkdir /opt/mcr && \
    wget --quiet -P /opt/mcr_install http://ssd.mathworks.com/supportfiles/downloads/R2020a/Release/0/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2020a_glnxa64.zip && \
    unzip -q /opt/mcr_install/MATLAB_Runtime_R2020a_glnxa64.zip -d /opt/mcr_install && \
    /opt/mcr_install/install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    rm -rf /opt/mcr_install /tmp/*

# Configure environment
ENV MCR_VERSION v98
ENV LD_LIBRARY_PATH /usr/lib/fsl/5.0:/opt/mcr/${MCR_VERSION}/runtime/glnxa64:/opt/mcr/${MCR_VERSION}/bin/glnxa64:/opt/mcr/${MCR_VERSION}/sys/os/glnxa64:/opt/mcr/${MCR_VERSION}/sys/opengl/lib/glnxa64
ENV MCR_INHIBIT_CTF_LOCK 1
ENV MCRPath /opt/mcr/${MCR_VERSION}

# Configure DPABI
RUN mkdir /opt/ENIGMA_RMD_Vertex_Stats_StandAlone
COPY ./ENIGMA_RMD_Vertex_Stats_StandAlone/* /opt/ENIGMA_RMD_Vertex_Stats_StandAlone/
COPY ./license.txt /opt/freesurfer/license.txt
RUN chmod 777 -R /opt/ENIGMA_RMD_Vertex_Stats_StandAlone/

# Extract ctf for singularity support
RUN /opt/ENIGMA_RMD_Vertex_Stats_StandAlone/run_ENIGMA_RMD_Vertex_Stats_StandAloneFixModel7.sh /opt/mcr/${MCR_VERSION} || true

RUN chmod 777 -R /opt/ENIGMA_RMD_Vertex_Stats_StandAlone/


# For freeview
RUN apt-get install -y qt4-default libjpeg62


ENTRYPOINT ["/opt/ENIGMA_RMD_Vertex_Stats_StandAlone/run_ENIGMA_RMD_Vertex_Stats_StandAloneFixModel7.sh", "/opt/mcr/v98"]


