foreach $t (580,600,620,640,660,680) {
	foreach $q ("a","b") {
		open A,"cm$t$q.fin" or die;
		open B,">temp";
		while(<A>) {
			s/19000/21000/ if /TZ/;
			s/850/900/ if /chi_inf/;
			print B;
		}
		close A;
		close B;
		`mv temp cm$t$q.fin`;
	}
}
