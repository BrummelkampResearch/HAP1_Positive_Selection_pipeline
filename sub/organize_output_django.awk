BEGIN {
	FS="\t"
}
{
	Gene=$1
	CT=$2
	TCT=$3
	NM=$4
	TNM=$5
	Pv=$6
	FCPv=$7
	RADIUS=$8
	SEQ=$9
	if (CT == 0) {
		CCT=1;
		CCTC=TCT-1;
	}
	else {
		CCT=CT;
		CTCT=TCT;
	}
	TI=NM+CCT
	MI=((NM/TNM)/(CCT/CTCT))
	print ","Gene","CT","TCT","CCT","CTCT","NM","TNM","Pv","FCPv","TI","MI","RADIUS","SEQ
}