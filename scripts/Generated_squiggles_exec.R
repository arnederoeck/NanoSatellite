args=commandArgs(trailingOnly = TRUE)

squiggle_file=args[1]
prefix=args[2]
path=args[3]
current_dir=args[4]

prefix_split=strsplit(prefix,"/")
set_title=prefix_split[[1]][length(prefix_split[[1]])]

report_path=paste(path,"scripts/Generated_squiggles_graph_report.Rmd",sep="/")

rmarkdown::render(report_path,params=list(squiggle_file=squiggle_file,set_title=set_title,current_dir=current_dir),output_file = paste(current_dir,"/",prefix,"_theoretical_squiggles.html",sep=""),quiet=TRUE)
