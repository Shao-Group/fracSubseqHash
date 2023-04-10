#!/usr/bin/env k8

/*
  Given a file of efa headers (see grepMappedReads.js) sorted by -k6,6 -k7,7n -k8,8n, output all pairs whose
  mappings overlap at least min_ovlp on the reference genome.

  Based on paftools.js in minimap2.

  By: Ke@PSU
  Last edited: 10/19/2022
*/


/*******************************
 * Command line option parsing *
 * from paftools.js
 *******************************/

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}


function main(args)
{
    var c, min_ovlp = 15;
    while ((c = getopt(args, "l:")) != null) {
	if (c == 'l') min_ovlp = parseInt(getopt.arg);
    }
    if (args.length - getopt.ind < 1) {
	print("Usage: genGroundTruth.js [options] <headers.efa>");
	print("Options:");
	print("  -l INT     min overlap length [15]");
	exit(1);
    }

       
    var buf = new Bytes();
    var fin = new File(args[getopt.ind]);
    var fout = new File(args[getopt.ind]+".pairs", "w");

    var a=[];
    while (fin.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	var name = parseInt(t[0].substring(1));
	var ctg = t[5], st = parseInt(t[6]), ed = parseInt(t[7]);
	while (a.length > 0){
	    if(a[0][0] == ctg && a[0][2] > st) break;
	    else a.shift();
	}
	for (var i=0; i<a.length; ++i){
	    var len = (ed > a[i][2] ? a[i][2] : ed) - st;
	    if (len >= min_ovlp){
		fout.write((a[i][3] < name ? a[i][3]+" "+name : name+" "+a[i][3])+" "+len+"\n");
	    }
	}
	a.push([ctg, st, ed, name]);
    }
    fin.close();
    fout.close();
    buf.destroy();
}

main(arguments);
