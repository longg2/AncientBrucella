#!/usr/bin/awk -f
BEGIN{
	FS="\t"
	OFS="\t"
	#print "File", "Fragment", "Mean", "Std. Dev", "CV","Percent Coverage"
	}

/group_/{
	print $7
	}
!/group_/{
	print $1
}
