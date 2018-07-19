library(plotrix)
draw3Dcircle = function(x,y,r,hue,sat=.5,num=30) {
	if (num==0) {
		if (hue == "grey") {
			col=hsv(.6,.05,.8)
		} else {
			col=hsv(hue,sat,1)
		}
		draw.circle(x,y,r,col=col,border=FALSE)
	} else {
	for (i in 0:num) {
		if (sat<i/((110/30)*num)) {
			sat2 = sat
		} else {
			sat2 = i/((110/30)*num)
		}
		if (hue == "grey") {
			col=hsv(.6,.1,.7+i/((200/30)*num))
		} else {
			col=hsv(hue,sat-sat2,1)
		}
		draw.circle(x+(1/(sqrt(2)))*r*(i/(num*2)),y+(1/(sqrt(2)))*r*(i/(num*2)),r*(1-i/(num*4/3)),col=col,border=FALSE)
	}
	}

}
