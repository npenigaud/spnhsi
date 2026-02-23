export ROOTDIR=${HOME}/install-newfxtranacdc

PATH="$ROOTDIR/perl5/bin${PATH:+:${PATH}}"; export PATH;
PERL5LIB="$ROOTDIR/perl5/lib/perl5${PERL5LIB:+:${PERL5LIB}}"; export PERL5LIB;
PERL_LOCAL_LIB_ROOT="$ROOTDIR/perl5${PERL_LOCAL_LIB_ROOT:+:${PERL_LOCAL_LIB_ROOT}}"; export PERL_LOCAL_LIB_ROOT;
PERL_MB_OPT="--install_base \"$ROOTDIR/perl5\""; export PERL_MB_OPT;
PERL_MM_OPT="INSTALL_BASE=$ROOTDIR/perl5"; export PERL_MM_OPT;

export PATH=${HOME}/install-newfxtranacdc/fxtran-acdc_full/bin:$PATH
export PATH=${HOME}/install-newfxtranacdc/fxtran/bin:$PATH


source /home/afar/modules/use.sh
module load rocm
module load afar/22.2.0





