run("Subtract Background...", "rolling=25");
run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");
setOption("BlackBackground", false);
/*
run("Make Binary");
run("Invert");
*/
run("Grid ", "grid=Lines area=10000 color=Cyan");
