#####################################################
#         Olivier Catoni & Gautier Appert           #
#         GNU LESSER GENERAL PUBLIC LICENSE         # 
#####################################################


# require(tiff)
# require(jpeg)
# require(gtools)
# library(gWidgets2)
# options(guiToolkit="RGtk2")
# library(Rcpp)

paint <- function(W) { 
    WsIndex = sort(unique(W+1)) 
    number_of_colors = length(WsIndex)
    WsColors = rep(0, max(WsIndex))
    WsColors[WsIndex] = 1 : number_of_colors
    if (number_of_colors > 1) {
      colorIndex = cbind(matrix(0,3,1), col2rgb(rainbow(number_of_colors-1)) / 255)
      return(t(colorIndex[,WsColors[W+1]]))
    } else {
      return(matrix(0,length(W), 3))
    }
}

showPatchImage <- function(h, ...) {

  dev.set(3) # the training set and results display window

  if (firstNewLabel == 0) {# returns if no patches are computed yet 
    return()
  }

	if (is.na(svalue(patch_image_nb))) { 
		return()
	}

	if (svalue(syntax_level) < 2) { 
		return()
	} else if (svalue(syntax_level) == 2) { 
    dim(Wsx) <<- dim(x)
		image_patches = sort(unique(Wsx[,svalue(patch_image_nb)]))+1
		if (length(image_patches) > 1) { 
    	patch_nb[] <<- image_patches 
		} else { 
    	patch_nb[] <<- c(image_patches, NA) 
		}
		mask = Wsx[,svalue(patch_image_nb)] != (svalue(patch_nb)-1)
		w2 = x[,svalue(patch_image_nb)] %o% c(1,1,1) 
		w2[mask,] = rep(1,sum(mask)) %o% background 
    dim(w2) = c(mpix$size, 3)
    # dim(w2) = mpix$size
  } else if ((svalue(syntax_level) %% 2) == 1) {
    current_s_level = (svalue(syntax_level)-1) %/% 2
		Twx <<- Twx_l[[current_s_level]]
    dim(Twx) <<- c(number_of_patches, dim(x)[2])
    dim(Wsx) <<- dim(x)
		image_patches = sort(unique(Twx[,svalue(patch_image_nb)]))+1
		if (image_patches[1] == 0) {
			image_patches = image_patches[-1]
		}
		if (length(image_patches) > 1) { 
    	patch_nb[] <<- image_patches 
		} else  {
    	patch_nb[] <<- c(image_patches, NA)
		}
		mask = (Twx[,svalue(patch_image_nb)] != (svalue(patch_nb)-1))[
			Wsx[,svalue(patch_image_nb)]+1]
		w2 = x[,svalue(patch_image_nb)] %o% c(1,1,1) 
		w2[mask,] = rep(1,sum(mask)) %o% background 
    dim(w2) = c(mpix$size, 3)
  } else if ((svalue(syntax_level) %% 2) == 0) {
    current_s_level = (svalue(syntax_level)-1) %/% 2 
		Uwx <<- Uwx_l[[current_s_level]]
    dim(Uwx) <<- c(number_of_patches, dim(x)[2])
    dim(Wsx) <<- dim(x)
		image_patches = sort(unique(Uwx[,svalue(patch_image_nb)]))+1
		if (image_patches[1] == 0) {
			image_patches = image_patches[-1]
		}
		if (length(image_patches) > 1) { 
    	patch_nb[] <<- image_patches 
		} else  {
    	patch_nb[] <<- c(image_patches, NA)
		}
		mask = (Uwx[,svalue(patch_image_nb)] != (svalue(patch_nb)-1))[
			Wsx[,svalue(patch_image_nb)]+1]
		w2 = x[,svalue(patch_image_nb)] %o% c(1,1,1) 
		w2[mask,] = rep(1,sum(mask)) %o% background 
    dim(w2) = c(mpix$size, 3)
  }

  w = x[, svalue(patch_image_nb)]
  dim(w) = mpix$size
  rasterImage(w,dim(w)[2],0,2*dim(w)[2],dim(w)[1], interpolate = FALSE)

	if (exists("w2")) { 
  	rasterImage(w2,0,0,dim(w2)[2],dim(w2)[1], interpolate = FALSE)
	}
}

showPatch <- function(h, ...) {

  dev.set(3) # the training set and results display window

  if (firstNewLabel == 0) {# returns if no patches are computed yet 
    return()
  }

	if (is.na(svalue(patch_nb))) {
		return()
	}

	if (svalue(syntax_level) < 2) { 
		return()
	} else if (svalue(syntax_level) == 2) { 
    dim(Wsx) <<- dim(x)
		patch_images = which(colSums(Wsx == (svalue(patch_nb)-1)) > 0)
		if (length(patch_images) > 1) { 
    	patch_image_nb[] <<- patch_images 
		} else  { 
    	patch_image_nb[] <<- c(patch_images, NA) 
		}
		w2 = x[,svalue(patch_image_nb)] %o% c(1,1,1) 
		mask = Wsx[,svalue(patch_image_nb)] != (svalue(patch_nb)-1)
		w2[mask,] = rep(1,sum(mask)) %o% background 
    dim(w2) = c(mpix$size, 3)
  } else if ((svalue(syntax_level) %% 2) == 1) {
    current_s_level = (svalue(syntax_level)-1) %/% 2
		Twx <<- Twx_l[[current_s_level]]
    dim(Twx) <<- c(number_of_patches, dim(x)[2]) 
    dim(Wsx) <<- dim(x)
		patch_images = which(colSums(Twx == (svalue(patch_nb)-1)) > 0)
		if (length(patch_images) > 1) { 
    	patch_image_nb[] <<- patch_images 
		} else if (length(patch_images) <= 1) {
			patch_image_nb[] <<- c(patch_images, NA)
		}
		w2 = x[,svalue(patch_image_nb)] %o% c(1,1,1) 
		mask = (Twx[,svalue(patch_image_nb)] 
			!= (svalue(patch_nb)-1))[Wsx[,svalue(patch_image_nb)]+1]
		if (sum(mask) > 0) { 
			w2[mask,] = rep(1,sum(mask)) %o% background 
		}
    dim(w2) = c(mpix$size, 3)
  } else if ((svalue(syntax_level) %% 2) == 0) {
    current_s_level = (svalue(syntax_level)-1) %/% 2
		Uwx <<- Uwx_l[[current_s_level]]
    dim(Uwx) <<- c(number_of_patches, dim(x)[2]) 
    dim(Wsx) <<- dim(x)
		patch_images = which(colSums(Uwx == (svalue(patch_nb)-1)) > 0)
		if (length(patch_images) > 1) { 
    	patch_image_nb[] <<- patch_images 
		} else if (length(patch_images) <= 1) {
			patch_image_nb[] <<- c(patch_images, NA)
		}
		w2 = x[,svalue(patch_image_nb)] %o% c(1,1,1) 
		mask = (Uwx[,svalue(patch_image_nb)] 
			!= (svalue(patch_nb)-1))[Wsx[,svalue(patch_image_nb)]+1]
		if (sum(mask) > 0) { 
			w2[mask,] = rep(1,sum(mask)) %o% background 
		}
    dim(w2) = c(mpix$size, 3)
  }

  w = x[,svalue(patch_image_nb)]
  dim(w) = mpix$size
  rasterImage(w,dim(w)[2],0,2*dim(w)[2],dim(w)[1], interpolate = FALSE)

	if (exists("w2")) { 
  	rasterImage(w2,0,0,dim(w2)[2],dim(w2)[1], interpolate = FALSE)
	}
}


# Display of training set samples with the computed patches
showImage <- function(h, ...) {

	dev.set(3) # the training set and results display window
  w = x[,svalue(imagenb)]
  dim(w) = mpix$size
	if (firstNewLabel > 0) { 
  	rasterImage(w,dim(w)[2],0,2*dim(w)[2],dim(w)[1], interpolate = FALSE)
	} else {# returns if no patches are computed yet 
  	rasterImage(w,0,0,dim(w)[2],dim(w)[1], interpolate = FALSE)
		return()
	}


	if (svalue(syntax_level) == 1) {
		dim(Wsx) <<- dim(lx)
		w2 = exp(patch_mean[cbind(1:dim(lx)[1],Wsx[,svalue(imagenb)]+1)]) - epsilon
		dim(w2) = mpix$size
	} else if (svalue(syntax_level) == 0) {
  	w2 = x[,svalue(imagenb)]
  	dim(w2) = mpix$size
	} else if (svalue(syntax_level) == 2) { 
		dim(Wsx) <<- dim(lx)
		image_patches = sort(unique(Wsx[,svalue(imagenb)]))+1
		if (length(image_patches) > 1) { 
			patch_nb[] <<- image_patches
		} else {
			patch_nb[] <<- c( image_patches, NA)
		}
		patch_images = which(colSums(Wsx == (svalue(patch_nb)-1)) > 0)
		if (length(patch_images) > 1) { 
    	patch_image_nb[] <<- patch_images 
		} else { 
    	patch_image_nb[] <<- c( patch_images, NA) 
		}
    dim(Wsx) <<- dim(x)
    w2 = paint(Wsx[,svalue(imagenb)])
    dim(w2) = c(mpix$size, 3)
  } else if ((svalue(syntax_level)%%2) == 1) {
    current_s_level = (svalue(syntax_level)-1) %/% 2 
    dim(Wsx) <<- dim(x)
		Twx <<- Twx_l[[current_s_level]]
		dim(Twx) <<- c(number_of_patches, dim(x)[2])
		image_patches = sort(unique(Twx[,svalue(imagenb)]))+1
		if (image_patches[1] == 0) {
			image_patches = image_patches[-1]
		}
		if (length(image_patches) > 1) { 
			patch_nb[] <<- image_patches 
		} else {
			patch_nb[] <<- c(image_patches, NA)
		}
		patch_images = which(colSums(Twx == (svalue(patch_nb)-1)) > 0)
		if (length(patch_images) > 1) { 
			patch_image_nb[] <<- patch_images 
		} else if (length(patch_images) == 1) {
			patch_image_nb[] <<- c(patch_images,NA)
		} 
    w2 = paint(Twx[Wsx[,svalue(imagenb)]+1,svalue(imagenb)])
    dim(w2) = c(mpix$size, 3)
  } else if ((svalue(syntax_level) %% 2) == 0) {
    current_s_level = (svalue(syntax_level) - 1) %/% 2
    dim(Wsx) <<- dim(x)
		Uwx <<- Uwx_l[[current_s_level]]
		dim(Uwx) <<- c(number_of_patches, dim(x)[2])
		image_patches = sort(unique(Uwx[,svalue(imagenb)]))+1
		if (image_patches[1] == 0) {
			image_patches = image_patches[-1]
		}
		if (length(image_patches) > 1) { 
			patch_nb[] <<- image_patches 
		} else {
			patch_nb[] <<- c(image_patches, NA)
		}
		patch_images = which(colSums(Uwx == (svalue(patch_nb)-1)) > 0)
		if (length(patch_images) > 1) { 
			patch_image_nb[] <<- patch_images 
		} else {
			patch_image_nb[] <<- c(patch_images,NA)
		}
    w2 = paint(Uwx[Wsx[,svalue(imagenb)]+1,svalue(imagenb)])
    dim(w2) = c(mpix$size, 3)
  }
 
	if (exists("w2")) {
  	rasterImage(w2,0,0,dim(w2)[2],dim(w2)[1], interpolate = FALSE) 
	}
}

# resets myBeta parameter in conjunction with computePatch 
resetScale <- function(h, ...) { 
	myBeta <<- 2*10^svalue(scaleParameter) # updates the value of myBeta
} 

quitFunction <- function(h, ...) { 
if (FALSE) { 
	if (firstNewLabel != 0) { 
		save(patch_mean, patch_variance, patch_weights, Crit, Asets, Bsets, myBeta, myThreshold, kvalue, 
			firstNewLabel, lastNewLabel, firstMerge, file = "dataDump.RData") 
	} else { # firstNewLabel == 0
		save( firstNewLabel, firstMerge, file = "dataDump.RData") 
	}
}
	quit(save = "no")
} 

patchNumber1 <- function(h, ...) { # called when the number of patches to be 
	# computed is changed 
	kvalue[1] <<- strtoi(svalue(h$obj)) 
}

patchNumber2 <- function(h, ...) { # called when the number of patches to be 
	# computed is changed 
	kvalue[2] <<- strtoi(svalue(h$obj)) 
}

patchNumber3 <- function(h, ...) { # called when the number of patches to be 
	# computed is changed 
	kvalue[3] <<- strtoi(svalue(h$obj)) 
}

createSample <- function(h, ...) { 

sampleSize <<- svalue(sampleSlider) 
positions = c(
	sample((mpix$frame[1] - mpix$start[1] - mpix$size[1])%/% mpix$steps[1], 
		sampleSize, replace = TRUE),
	sample((mpix$frame[2] - mpix$start[2] - mpix$size[2])%/% mpix$steps[2], 
		sampleSize, replace = TRUE) 
)
dim(positions) = c(sampleSize, 2)

if (exists("x")) { 
	x <<- c(x, 
	sapply(1:sampleSize, function(k) {
	return(
		newImage[mpix$start[1] + mpix$steps[1] * positions[k,1] + (1:mpix$size[1]), 
			mpix$start[2] + mpix$steps[1] * positions[k,2] + (1:mpix$size[2])]
	)}
	))
} else {
	x <<-  
	sapply(1:sampleSize, function(k) {
	return(
		newImage[mpix$start[1] + mpix$steps[1] * positions[k,1] + (1:mpix$size[1]), 
			mpix$start[2] + mpix$steps[1] * positions[k,2] + (1:mpix$size[2])]
	)}
	)
}
	
dim(x) <<- c(prod(mpix$size),length(x)/prod(mpix$size))

x <<- x[,sample(dim(x)[2])]

imagenb[] <<- seq(1, dim(x)[2], by=1)

cat("dim(x) = ", dim(x), "\n")

epsilon <<- 10^(svalue(epsilonSlider))
lx <<- log ( x + epsilon ) 

# End of creation of the image sample x. 

dev.off(3)
# Display training sample and / or results
x11(title="training set",xpos = -1, ypos=200, width=6, height=6)
par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
plot(c(0, mpix$size[2]), c(0, mpix$size[1]), type = "n",
  ann = FALSE, axes = FALSE, asp = 1)

	firstNewLabel <<- 0 # initialization of patch numbers
} 

computePatches <- function(h, ...) { 
	# computes patches number firstNewLabel 
	# to lastNewLabel
	
	if (!firstMerge) {
		cat("You cannot split again once you merged patches !\n")
		return()
	}

	s_level <<- 0
		firstNewLabel <<- dim(lx)[2] + 1
		lastNewLabel <<- dim(lx)[2] + kvalue[1] 
		patch_mean <<- cbind(lx, matrix(0, dim(lx)[1], kvalue[1]))
		cat("dim(patch_mean) = ", dim(patch_mean), "\n")
		patch_variance <<- matrix(0, dim(lx)[1], lastNewLabel)
		patch_weights <<- rep(1, lastNewLabel)
		Crit <<- rep(0, lastNewLabel)
		myBeta <<- 2*10^svalue(scaleParameter)
		Crit[1:dim(lx)[2]] <<- myBeta * dim(lx)[1] - 1 
		# the criterion is now myBeta*|B| - |A| 
		Asets <<- integer(dim(lx)[2]*lastNewLabel) 
		dim(Asets) <<- c(dim(lx)[2], lastNewLabel)
		Bsets <<- integer(dim(lx)[1] * lastNewLabel)
		dim(Bsets) <<- c(dim(lx)[1], lastNewLabel)
		Wsx <<- integer(length(lx))
		dim(Wsx) <<- dim(lx)
		myThreshold <<- 0

	cat("Computes patches number ", firstNewLabel, 
			" to ", lastNewLabel, "\n")
	number_of_patches <<- lastNewLabel

# loads the parameter list to be passed to C++
	paramList <<- list ( n = dim(lx)[2], d = dim(lx)[1], 
		betaCoeff = myBeta, # loaded above 
		thresholdCoeff = svalue(thresholdSlider)^2,
		thresholdValue = myThreshold,
		firstNewLabel = firstNewLabel,
		lastNewLabel = lastNewLabel
	) 
	cat("computeLoop !\n")
	computeLoop(patch_mean,patch_variance,patch_weights,Crit, Asets, Bsets, 
		Wsx, paramList) 
	myThreshold <<- paramList$thresholdValue # updates the value of 
	# myThreshold from the value computed by computeLoop
	insert(messages,paste("sqrt(threshold) = ", sqrt(myThreshold)))
	cat ("sqrt(threshold) = ", sqrt(myThreshold), "\n")

	syntax_level[] <<- seq(0,2,by=1)
	Twx_l <<- list()
	Uwx_l <<- list()
	number_of_patches_l <<- list(number_of_patches)
	number_of_syntax_labels_l <<- list() 

	dev.off(4)
	# Display of patch labels 
	x11(title="Patch labels", width=16, 
		# height=16*dim(Asets)[1]/(dim(Asets)[2]-dim(Asets)[1]))
		height=8)
	par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
	plot(c(0,16), c(0,8), type = "n", 
			ann = FALSE, axes = FALSE)
	rasterImage(Asets[,order(colSums(Asets),decreasing=TRUE)],0,0, 16, 8, interpolate = FALSE)

	dev.off(3)
	# Display training sample and / or results
	x11(title="training set",xpos = -1, ypos=200, width=12, height=6)
	par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
	plot(c(0, 2*mpix$size[2]), c(0, mpix$size[1]), type = "n",
  	ann = FALSE, axes = FALSE, asp = 1)

}

iteratePatches <- function(h, ...) { 
	if (firstNewLabel == 0) { 
		cat("You should create patches first.\n")
		return()
	}
	s_level <<- s_level + 1 # s_level == 1 on first call
	number_of_patches_l[[s_level]] <<- ifelse(s_level == 1, number_of_patches, 
		number_of_syntax_labels_l[[s_level-1]])
	current_number_of_patches <<- number_of_patches_l[[s_level]]

	firstNewLabel <<- current_number_of_patches + 1
  lastNewLabel <<- current_number_of_patches + kvalue[2]
	number_of_merged_labels <<- lastNewLabel
	cat("First merge !\n")
	Asets <<- c(Asets, integer(dim(lx)[2]*kvalue[2]))
	dim(Asets) = c(dim(lx)[2], number_of_merged_labels)
	cat("dim(Asets) = ", dim(Asets), "\n")
	Asizes <<- integer(lastNewLabel)
	Bsets <<- integer( current_number_of_patches * number_of_merged_labels) 
	Bsizes <<- integer(number_of_merged_labels)
	dim(Bsets) <<- c(current_number_of_patches, number_of_merged_labels)
	Twx <<- integer(current_number_of_patches*dim(lx)[2])
	dim(Twx) <<- c(current_number_of_patches, dim(lx)[2])
	Jt <<- integer(2*number_of_merged_labels)
	dim(Jt) <<- c(2, number_of_merged_labels)  
	paramList <<- list ( n = dim(lx)[2], d = dim(lx)[1], # d is not necessary ? 
		betaCoeff = myBeta, 
		firstNewLabel = firstNewLabel,
		lastNewLabel = lastNewLabel,
		firstMerge = firstMerge
	) 
	cat("computeMerge !\n")
	
	computeMerge(Asets, Asizes, Bsets, Bsizes, Twx, Jt, paramList)

	firstNewLabel <<- lastNewLabel+1
	lastNewLabel <<- lastNewLabel + kvalue[3]
	number_of_syntax_labels <<- lastNewLabel

	Cont <<- integer(number_of_merged_labels*number_of_syntax_labels)
	ContOrder <<- integer(number_of_merged_labels)
	dim(Cont) <<- c(number_of_merged_labels,number_of_syntax_labels)
	computeContext(Jt, Cont, ContOrder, current_number_of_patches, number_of_merged_labels)

	# Display of syntax labels 
	x11(title="Context matrix", width=8, height=8)
	par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
	plot(c(0,12), c(0,12), type = "n", 
			ann = FALSE, axes = FALSE)
	dim(Cont) = c(number_of_merged_labels, number_of_syntax_labels)
	rasterImage(Cont[,1:number_of_merged_labels],0,0, 12, 12, interpolate = FALSE)

	Contsizes <<- integer(number_of_syntax_labels)
	Bsets <<- integer( number_of_merged_labels * number_of_syntax_labels) 
	Bsizes <<- integer(number_of_syntax_labels)
	dim(Bsets) <<- c(number_of_merged_labels, number_of_syntax_labels)
	Gtt <<- integer(number_of_merged_labels^2)
	dim(Gtt) <<- c(number_of_merged_labels, number_of_merged_labels)
	Jg <<- integer(2*number_of_syntax_labels)
	dim(Jg) <<- c(2, number_of_syntax_labels)  
	paramList <<- list ( 
		n = number_of_merged_labels, 
		d = number_of_merged_labels, 
		betaCoeff = myBeta, 
		firstNewLabel = firstNewLabel,
		lastNewLabel = lastNewLabel,
		firstMerge = firstMerge # not used any more
	) 
	
	computeMerge(Cont, Contsizes, Bsets, Bsizes, Gtt, Jg, paramList)

if (FALSE) { 
	# Display of syntax labels 
	x11(title="Context matrix", width=8, height=8)
	par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
	plot(c(0,12), c(0,12), type = "n", 
			ann = FALSE, axes = FALSE)
	dim(Cont) = c(number_of_merged_labels, number_of_syntax_labels)
	rasterImage(Cont,0,0, 12, 12, interpolate = FALSE)
}

	Gt <<- integer(number_of_merged_labels)
	Uwx <<- integer(current_number_of_patches * dim(lx)[2])
	Axu <<- integer(dim(lx)[2] * lastNewLabel)
	applySyntax(Twx, Gtt, Jt, Gt, Uwx, Asets, Axu, 
		number_of_syntax_labels, 
		number_of_merged_labels, 
		current_number_of_patches, 
		dim(lx)[2]) 

cat("applySyntax done !\n")

	firstMerge <<- FALSE
	cat("pairs = \n")
	print(t(Jt[,(current_number_of_patches+1):number_of_merged_labels])+1)

# syntax analysis 
	syntax_level[] <<- seq(from = 0, to = s_level*2 + 2, by = 1) 
	number_of_syntax_labels_l[[s_level]] <<- number_of_syntax_labels
	if (s_level > 1) { 
		Twx_l[[s_level]] <<- integer(number_of_patches * dim(lx)[2])
		combineLabels(Twx_l[[s_level]], Twx, Uwx_l[[s_level-1]],
			number_of_patches, current_number_of_patches, dim(lx)[2])
		Uwx_l[[s_level]] <<- integer(number_of_patches * dim(lx)[2])
		combineLabels(Uwx_l[[s_level]], Uwx, Uwx_l[[s_level-1]],
			number_of_patches, current_number_of_patches, dim(lx)[2])
	} else {
		Twx_l[[s_level]] <<- Twx
		Uwx_l[[s_level]] <<- Uwx
	}
	Asets <<- Axu
}

showPatchFunction <- function(h, ...) {
	patch_nb$invoke_change_handler()
}

showAllFunction <- function(h, ...) {
	if (svalue(imagenb) != svalue(patch_image_nb)) {
		svalue(imagenb) <<- svalue(patch_image_nb)
	} else {
		imagenb$invoke_change_handler()
	}
}

showFrame <- function(h, ...) { 

	dev.set(2)
	mpix$start[1] <<- svalue(ULFrameX)
	mpix$start[2] <<- svalue(ULFrameY)
	mpix$frame[1] <<- svalue(LRFrameX)
	mpix$frame[2] <<- svalue(LRFrameY)
	w = newImage %o% c(1,1,1)
	w[mpix$start[1]:mpix$frame[1], mpix$start[2]:mpix$frame[2],] = 
	w[mpix$start[1]:mpix$frame[1], mpix$start[2]:mpix$frame[2],1] %o% c(1,1,0)
	rasterImage(w,0,0,dim(w)[2],dim(w)[1], interpolate = FALSE)
}

loadImage <- function(h, ...) {

	fileName = gfile("Image file name", type = "open")
	newImage <<- readJPEG(fileName)[,,2] 

	ULFrameX[] <<- seq(1, dim(newImage)[1], by=1)
	svalue(ULFrameX) <<- floor(dim(newImage)[1] * 0.2)
	ULFrameY[] <<- seq(1, dim(newImage)[2], by=1)
	svalue(ULFrameY) <<- floor(dim(newImage)[2] * 0.2)
	LRFrameX[] <<- seq(1, dim(newImage)[1], by=1)
	svalue(LRFrameX) <<- floor(dim(newImage)[1] * 0.8)
	LRFrameY[] <<- seq(1, dim(newImage)[2], by=1)
	svalue(LRFrameY) <<- floor(dim(newImage)[2] * 0.8)

	# Display of original images 
	cat("dev number = ", dev.cur(), "\n")
	dev.off(2)
	x11(title="Original images", width=8, 
		height=8*dim(newImage)[1]/dim(newImage)[2])
	par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
	plot(c(0,dim(newImage)[2]), c(0,dim(newImage)[1]), type = "n",
  	ann = FALSE, axes = FALSE, asp = 1)
	rasterImage(newImage,0,0,dim(newImage)[2],dim(newImage)[1], 
		interpolate = FALSE)
}

PatchProcessWidget <- function() { 

kvalue <<- rep(0,4)

win <<- gwindow("Prepare sample") # command window and main title
gp <<- ggroup(horizontal=FALSE, cont=win) # container

tmp <<- gbutton("Load new image", cont=gp, handler=loadImage,
expand=TRUE)

tmp <<- ggroup(cont = gp, expand = TRUE)
tmp2 <<- glabel(" Upper left frame corner : ",cont=tmp)
ULFrameX <<- gslider(from=1, to=2,
by=1, value=150, cont=tmp, expand=TRUE) 
ULFrameY <<- gslider(from=1, to=2,
by=1, value=350, cont=tmp, expand=TRUE) 
tmp <<- ggroup(cont = gp, expand = TRUE)
tmp2 <<- glabel(" Lower right frame corner :",cont=tmp)
LRFrameX <<- gslider(from=1, to=2,
by=1, value=800, cont=tmp, expand=TRUE) 
LRFrameY <<- gslider(from=1, to=2,
by=1, value=1200, cont=tmp, expand=TRUE) 
showFrameButton <<- gbutton("Set frame",cont=gp,handler=showFrame,
expand=TRUE)

tmp <<- gframe("Image size", cont=gp, expand=TRUE)
sizeSlider <<- gslider(from=1, to = 400, by=1, value = 300, cont=tmp, 
expand=TRUE)

tmp <<- ggroup(cont = gp, expand = TRUE)
tmp2 <<- gframe("Step size", cont=tmp, expand=TRUE)
stepSlider <<- gslider(from=1, to = 100, by=1, value = 5, cont=tmp2, 
expand=TRUE)
tmp2 <<- gframe("Sample size", cont=gp, expand=TRUE)
sampleSlider <<- gslider(from=1, to = 10000, by=1, value = 100, cont=tmp2, 
expand=TRUE)

# sets the image size
tmp <<- ggroup(cont = gp, expand = TRUE)

tmp2 <<- gframe("Epsilon value", cont=tmp, expand=TRUE)
epsilonSlider <<- gslider(from=-4, to = 2, by=0.1, value = -2, cont=tmp2, expand=TRUE)

sizeButton <<- gbutton("Add to training set", handler = createSample, 
		cont = tmp, expand = TRUE)


win2 <<- gwindow("Compute syntax", parent= c(0,0)) # commands window and main title
gp <<- ggroup(horizontal=FALSE, cont=win2) # container

tmp <<- ggroup(cont = gp, expand = TRUE)

tmp2 <<- gframe("Beta", cont=tmp, expand=TRUE) # reads myBeta 
scaleParameter <<- gslider(from=-6, to = 4, by=0.1, value = 0, cont=tmp2,
handler = resetScale, expand=TRUE)
tmp2 <<- gframe("Threshold value", cont=tmp, expand=TRUE)
thresholdSlider <<- gslider(from=0, to = 2, by=0.01, value = 0.1, cont=tmp2, expand=TRUE)

# Sets the number of patches to be computed when pressing the Compute patches
# button 
tmp <<- ggroup(cont = gp, expand = TRUE)
tmp2 <<- glabel(" Number of patches : ", cont=tmp)
patchNumberFrame1 <<- gedit("1000", width = 4, handler = patchNumber1, cont=tmp)
addHandlerBlur(patchNumberFrame1, handler = patchNumber1)
patchNumberFrame1$invoke_change_handler()


# button to compute initial patches 
computePatchButton <<- gbutton("Compute initial patches", 
handler = computePatches, cont = gp, expand = TRUE)

tmp <<- ggroup(cont = gp, expand = TRUE)
tmp2 <<- glabel(" Number of merged labels: ", cont=tmp)
patchNumberFrame <<- gedit("300", width = 4, handler = patchNumber2, cont=tmp)
addHandlerBlur(patchNumberFrame, handler = patchNumber2)
patchNumberFrame$invoke_change_handler()

tmp <<- ggroup(cont = gp, expand = TRUE)
tmp2 <<- glabel(" Number of syntax labels: ", cont=tmp)
patchNumberFrame <<- gedit("300", width = 4, handler = patchNumber3, cont=tmp)
addHandlerBlur(patchNumberFrame, handler = patchNumber3)
patchNumberFrame$invoke_change_handler()

tmp <<- ggroup(cont = gp, expand = TRUE)
# button to merge patches
iteratePatchButton <<- gbutton("Compute next level patches", handler = iteratePatches, cont
= tmp, expand = TRUE)
tmp <<- gframe("Syntax level", cont=gp, expand=TRUE) 
syntax_level <<- gslider(from=0,to=1,by=1,value=0,cont=tmp,expand=TRUE,
  handler = showImage)

# slider through training sample indices, shows the patches present in a given 
# training sample
tmp <<- gframe("Show patches in image number", cont=gp, expand=TRUE) 
imagenb <<- gslider(from=1,to=100,by=1,value=1,cont=tmp,expand=TRUE,
  handler = showImage)

tmp2 <<- ggroup(cont = gp, expand = TRUE)
tmp <<- gframe("Show one patch", cont=tmp2, expand=TRUE)
patch_nb <<- gslider(from=c(1,2),value=1,cont=tmp,expand=TRUE,
  handler = showPatch)
showAllButton <<- gbutton("Show all", handler = showAllFunction, 
	cont = tmp2, expand = FALSE)

tmp2 <<- ggroup(cont = gp, expand = TRUE)
tmp <<- gframe("Change image", cont=tmp2, expand=TRUE)
patch_image_nb <<- gslider(from=c(1,2), value=1,cont=tmp,expand=TRUE,
  handler = showPatchImage)
showPatchButton <<- gbutton("Show patch", handler = showPatchFunction, 
	cont = tmp2, expand = FALSE)

# quits the application
quitButton <<- gbutton("Quit", handler = quitFunction, cont = gp, expand = TRUE)

messages <<- gtext(cont = gp)

mpix <<- list ( 
  size = c(svalue(sizeSlider), svalue(sizeSlider)), 
			# c(300 , 300), # size of training windows 
  steps = c(svalue(stepSlider),svalue(stepSlider)), # vertical and horizontal steps taken to sample
		# training windows
  start = c(150, 350), # position of the first training window
  frame = c(800, 1200) # frame lower right corner
		# of training windows extraction area
) 

s_level <<- 0

background <<- c(0.05,0.4,0.05)

firstNewLabel <<- 0
firstMerge <<- TRUE

Sys.sleep(Inf) # to avoid quitting immediately
}


