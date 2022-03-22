collapseCI = function(pairs,SNP.gene) {
	count=0
	CIs = unique(c(pairs))
	CIs = CIs[CIs%in%SNP.gene[,2]==FALSE]
	done = c()
	for (i in 1:length(CIs)) {
		CI = CIs[i]
		if (CI %in% done) {
			next
		}
		seeds = c(pairs[pairs[,1]==CI|pairs[,2]==CI,])
		seeds = unique(seeds[seeds!=CI])
		matches=CI
		for (j in (i+1):length(CIs)) {
			if (j > length(CIs)) {
				break
			}
			CI2 = CIs[j]
			if (CI2 %in% done) {
				next
			}
			if (CI!=CI2) {
				seeds2 = c(pairs[pairs[,1]==CI2|pairs[,2]==CI2,])
				seeds2 = unique(seeds2[seeds2!=CI2])
				if (all(seeds2%in%seeds)&all(seeds%in%seeds2)) {
					matches = c(matches,CI2)
				}
			}
		}
		if (length(matches)>1) {
			count=count+1
			done = c(done,matches)
			matches = unique(matches)
			joint = paste(matches,collapse=",")
			for (match in matches) {
				pairs[pairs[,1]==match,1]=joint
				pairs[pairs[,2]==match,2]=joint
			}
		}
	}
	results=c()
	results$pairs = pairs
	results$count = count
	return(results)
}
