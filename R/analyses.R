#############
# RSTAN SPECS
#############
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###########
# FUNCTIONS
###########
niceThousands  <-  function (number) {
    formatC(number, format = 'fg', big.mark = ',')
}

fromCelsiusToKelvin  <-  function (temperatureInCelsius) {
	temperatureInCelsius + 273.15
}

readFile  <-  function (filePath, ...) {
	read.csv(filePath, header = TRUE, stringsAsFactors = FALSE, ...)
}

readWitrox  <-  function (filePath, ...) {
	read.table(filePath, skip = 12, header = TRUE, check.names = FALSE, sep = '\t', stringsAsFactors = FALSE, ...)
}

readAllWitroxes  <-  function (...) {	
	allWitroxes  <-  dir(full.names = TRUE, recursive = TRUE, ...)
	listOfFiles  <-  lapply(allWitroxes, readWitrox)
	names(listOfFiles)  <-  basename(allWitroxes)
	listOfFiles
}

cleanWitroxFiles  <-  function (witroxFiles) {
	lapply(witroxFiles, cleanRespirometryData)
}

cleanRespirometryData  <-  function (data) {
	dates  <-  data[['Date &Time [DD-MM-YYYY HH:MM:SS]']]
	dates  <-  lapply(dates, function (x) strsplit(x, ' ')[[1]])
	times  <-  sapply(dates, function (x) paste0(x[2:3], collapse = ' '))
	dates  <-  sapply(dates, function (x) x[1])

	oxSatCols  <-  grep('[% air saturation]', names(data), fixed = TRUE)
	newData    <-  cbind(date = dates, time = seq_along(times), data[, oxSatCols])
	toKeep     <-  sapply(newData, function (x) !all(is.na(x)))
	newData[, toKeep]
}

makeRespirationList  <-  function (cleanedWitroxFiles, verbose = TRUE) {
	out  <-  vector(mode = 'list', length = length(cleanedWitroxFiles))
	for (i in seq_along(out)) {
		message('\n\n', names(cleanedWitroxFiles)[i])
		out[[i]]  <-  analyseRespirometry(data = cleanedWitroxFiles[[i]], verbose = verbose, rescaleY = TRUE, hypoxiaThreshold = 0.8, alpha = 0.5, method = I('pc'))
	}
	names(out)  <-  names(cleanedWitroxFiles)
	out
}

analyseRespirometry  <-  function (data, verbose, rescaleY = TRUE, hypoxiaThreshold = 0.8, ...) {
	chaCols  <-  grep('CH ', names(data), fixed = TRUE)
	chaNums  <-  as.numeric(gsub('.+[CH]+[[:space:]]+([1-4])+.*', '\\1', names(data)[chaCols]))
	out         <-  vector(mode = 'list', length = length(chaNums))
	names(out)  <-  paste0('channel_', chaNums)
	for (i in seq_along(chaNums)) {
		if (verbose) {
			message(names(out)[i])
		}
		dat  <-  data[, c(2, chaCols[i])]
		if (rescaleY) {
			dat[, 2]  <-  dat[, 2] + (100 - dat[1, 2]) # rescale maximum value to 100% air saturation
		}
		dat  <-  dat[dat[, 2] >= hypoxiaThreshold, ] # exclude potential data that hits hypoxia threshold
		if (nrow(dat) > 150) {
			dat   <-  LoLinR::thinData(dat, xy = c(1, 2), by = round(nrow(dat) / 150))[[1]]
		}
		out[[i]]  <-  LoLinR::rankLocReg(xall = dat[, 1], yall = dat[, 2], ...)
	}
	out
}

cleanRawData  <-  function (rawData) {
	cols                 <-  c('20' = 'tomato', '26' = 'dodgerblue2', '29' = 'seagreen3', '32' = 'purple')
	shps                 <-  c('20' = 21, '26' = 22, '29' = 23, '32' = 24)
	rawData$tempCelsius  <-  as.numeric(sapply(rawData$tankID, function (z) strsplit(z, '_')[[1]][1]))
	rawData$tempKelvin   <-  fromCelsiusToKelvin(rawData$tempCelsius)
	rawData$boltzmannK   <-  8.62e-5
	rawData$invKT        <-  1 / rawData$boltzmannK * (1 / mean(rawData$tempKelvin, na.rm = TRUE) - 1 / rawData$tempKelvin)
	rawData$lnMass_g     <-  log(rawData$mass_g)
	rawData$color        <-  cols[as.character(rawData$tempCelsius)]
	rawData$shape        <-  shps[as.character(rawData$tempCelsius)]
	rawData$Date         <-  as.Date(gsub('/18', '/2018', rawData$Date, fixed = TRUE), format = '%d/%m/%Y')
	rawData$timeDays     <-  as.numeric(rawData$Date - as.Date('2018-04-12'))
	rawData$tankN        <-  sapply(rawData$tankID, function (z) strsplit(z, '_')[[1]][2])
	rawData
}

makeRespirationData  <-  function (respirationList, cleanedData) {
	data  <-  mapply(wrangleLoLinROutput, LoLinRChannelList = respirationList, fileName = names(respirationList), MoreArgs = list('cleanedData' = cleanedData, 'fullLoLinROutputList' = respirationList), SIMPLIFY = FALSE)
	data  <-  do.call(rbind.data.frame, data)
	data  <-  data[-grep('Control', data$tankID), ]
	data  <-  data[!is.na(data$respRate_mlO2_h), ]

	# convert metabolic rates in three steps:
	# first, from ml O2 / hour to mg O2 / h == (1.429 mg O2 / 1 ml O2) = 1 * 1.429
	# second, from mg O2 / h to g C / d == (24 h / 1 d) * (1 g / 1e3 mg) * (12 g C / 32 g O2) = 24 / 1e3 * 12 / 32 = 0.009
	# third, from g C / d to J / d == this conversion should be calculated based on the standard enthalpy combustion of glucose, i.e. a biochemical argument which corresponds to metabolism -- it includes both the energy to be sequestrated by ATP and the energy lost during this process. (2805 kJ / 1 mole Glucose) * (1 mole of Glucose / 6 moles C) * (1 mole C / 12 g C) = 2805/6/12 approx. 39 kJ per g C
	data$respRate_mgO2_h  <-  data$respRate_mlO2_h * 1.429
	data$respRate_gC_d    <-  data$respRate_mgO2_h  * 0.009
	data$respRate_j_d     <-  data$respRate_gC_d * 39e3
	data$lnRespRate_j_d   <-  log(data$respRate_j_d)
	rownames(data)  <-  NULL
	data
}

wrangleLoLinROutput  <-  function (LoLinRChannelList, fileName, cleanedData, fullLoLinROutputList) {
	# first select the right file from raw data and create channel labels
	# then find the right bacterial control, but if control, just ignore it
	# volume of chamber is 27 ml
	data                  <-  cleanedData[cleanedData$witroxFileName == fileName, ]
	data$ch               <-  channelLabelFromWitroxChannel(data$witroxChannel)
	control               <-  cleanedData[cleanedData$fishID == unique(data$controlChamber), ]
	controlCh             <-  channelLabelFromWitroxChannel(control$witroxChannel)
	controlResp           <-  fullLoLinROutputList[[control$witroxFileName]][[controlCh]]
	controlSlope          <-  median(controlResp$allRegs$b1)
	data$respRate_mlO2_h  <-  NA
	for (k in 1:nrow(data)) {
		subsetData   <-  LoLinRChannelList[[data$ch[k]]]$allRegs
		respSlope    <-  median(subsetData$b1)
		temperature  <-  strsplit(fileName, '_')[[1]][1]
		bo2          <-  c('20' = 6.4, '26' = 5.69, '29' = 5.38, '32' = 5.1)[temperature]
		data$respRate_mlO2_h[k]  <-  vo2(ma = respSlope, mb = controlSlope, v = 27 - data$mass_g[k], bo2 = bo2) # in ml O2 / hour
	}
	data
}

channelLabelFromWitroxChannel  <-  function (witroxChannel) {
	paste0('channel_', sapply(witroxChannel, substr, 1, 1))
}

vo2  <-  function (ma, mb = 0, v, bo2 = 5.80) {
	# ma is the rate of change (slope) of O2 saturation (in % / second)
	# mb is the rate of change (slope) of O2 saturation (in % / second) for the control, i.e. bacterial respiration. It is considered negligible by default, i.e. 0.
	# v is the volume of vial after removing the volume of the animal -- the density of fish weight here is considered to be 1 g of wet weight / 1 ml.
	# bo2 is the oxygen volume solubility (ml O2 / L) at a given temperature, default of 5.8 is for freshwater (salinity 0) at 25oC.
	# convert unit to ml O2 per hour
	bo2  <-  bo2 / 1e3 # convert to ml O2 / ml
	(-1 * ((ma - mb) / 100) * v * bo2) * 60 * 60 # in ml O2 / hour
}

cleanMitochondriaData  <-  function (fileName, data, respModel, growthModel) {
	###########
	# RATIONALE
	###########
	# The P/O ratio is derived from the O2 consumed (µmol / L) during state 3 respiration (i.e. between adding ADP to the solution until S4 starts). This quantity then divides the known amount of ADP added to the reaction (1 µL, original concentration of stock solution was 25 mM, i.e. 25e3 µmol / L)
	# The initial volume in the chamber is 220 µL, which is then filled up with 50 µL of actual tissue-derived mitochondria solution. Then, during S2, 5 µL of malate and 2.5 µL of piruvate are added before adding ADP at the beginning of S3, and then 1 µL of ADP is added; which means that total volume inside the chamber prior to S3 is: 220 µL + 50 µL + 5 µL + 2.5 µL + 1 µL = 278.5 µL

	###################
	# CALCULATION STEPS
	###################
	# 1 - determine O2 consumed (in µmol / L) from addition of ADP to chamber until start of S4 (i.e. S3)
	# 2 - convert O2 to O: multiply by 2
	# 3 - determine how much O is in final sample volume; e.g. final sample volume = 278.5 µl = 278.5e-6 L -> multiply [O] in µmol / L (from step 2 above) by 278.5e-6, which gives µmol of O in chamber (actual O used)
	# 4 - calculate how much ADP was added to the chamber. Original concentration of ADP was 25e3 µmol / L, and 1 µL (1e-6 L) was added. So, the total amount of ADP added to the chamber is 25e3 * 1e-6 = 0.025 µmol of ADP
	# 5 - divide the value from step 4 (here 0.025 µmol) by the O used (from 3) to get P/O ratio
	# 6 - The P added to the ADP is in the assay buffer and is assumed to be non-limiting, so the P in P/O ratio refers to one P added to ADP to make ATP, which is indicated by S3

	mito           <-  readFile(fileName)
	mito$consumed  <-  mito[['O2BeginningS3_O2_umol_l']]  - mito[['O2EndS3_O2_umol_l']]
	mito$poRatio   <-  0.025 / (mito$consumed * 2 * 278.5e-6)

	fixefs     <-  brms::fixef(respModel)
	ranR       <-  coef(respModel)$fishID
	ranG       <-  coef(growthModel)$fishID
	endData    <-  data[data$Date > as.Date('2018-07-13'), ] # keep only fish whose mass were measured one week after respirometry finished

	endData$alpha      <-  ranR[endData$fishID, 'Estimate', 'lnMass_g']
	endData$Er         <-  ranR[endData$fishID, 'Estimate', 'invKT']
	endData$boTs       <-  exp(ranR[endData$fishID, 'Estimate', 'Intercept'] + endData$Er * endData$invKT)
	endData$theta      <-  ranG[endData$fishID, 'Estimate', 'gAlpha_Intercept']
	endData$Es         <-  ranG[endData$fishID, 'Estimate', 'Es_Intercept']
	endData$Eg         <-  ranG[endData$fishID, 'Estimate', 'Eg_Intercept']
	endData$asympMass  <-  exp(ranG[endData$fishID, 'Estimate', 'lnAsympMass_Intercept'] + endData$Es * endData$invKT)
	endData$aoTs       <-  exp(ranG[endData$fishID, 'Estimate', 'lnGoTs_Intercept'] + endData$Eg * endData$invKT)

	endData$oxCon_j_d   <-  endData$boTs * endData$mass_g ^ endData$alpha
	# mR is estimated in J / d, back transform to mg O2 / h, and then to moles / day: 1 g of O2 has 0.031251171918947 moles
	endData$mR_molO2_d  <-  endData$oxCon_j_d / 39e3 / 0.009 * (0.031251171918947 * 1e-3) * 24 # in moles of O2 / day
	endData$poRatio     <-  mito$poRatio[match(endData$fishID, mito$fishID)]
	endData$molAtp_d    <-  endData$mR_molO2_d * endData$poRatio
	endData$growth_g_d  <-  endData$aoTs * endData$mass_g ^ endData$theta * (1 - (endData$mass_g / endData$asympMass) ^ (1 - endData$theta)) # estimate from growthModel
	endData$emAtp       <-  ((endData$molAtp_d * (1 - (endData$mass_g / endData$asympMass) ^ (1 - endData$alpha))) / endData$growth_g_d) * 1e3 # transform Em to mili moles / g
	endData$emOxC       <-  ((endData$oxCon_j_d * (1 - (endData$mass_g / endData$asympMass) ^ (1 - endData$alpha))) / endData$growth_g_d) / 1e3 # transform Em to kJ / g
	endData[!is.na(endData$molAtp_d), ] # 2 were excluded due to equipment mal-function
}

####################
# STATISTICAL MODELS
####################
putterGrowthModelBrms  <-  function (data) {
	data                <-  data[!is.na(data$mass_g) & !(data$fishID %in% c('32_2_5', '32_3_2')), ] # dead on week 2
	putterGrowthPriors  <-  c(brms::prior(normal(0, 5), nlpar = 'lnAsympMass'),
							brms::prior(normal(0, 5), nlpar = 'gAlpha'),
							brms::prior(normal(0, 5), nlpar = 'lnGoTs'),
	                        brms::prior(normal(0, 5), nlpar = 'Eg'),
							brms::prior(normal(0, 5), nlpar = 'Es'))

	putterGrowthInits    <-  list(list(lnAsympMass = log(0.4), gAlpha = 0.8, lnGoTs = log(1), Eg = -0.5, Es = 0.1),
	                            list(lnAsympMass = log(0.8), gAlpha = 0.75, lnGoTs = log(5), Eg = 0.8, Es = 1),
	                            list(lnAsympMass = log(1.0), gAlpha = 0.6, lnGoTs = log(0.1), Eg = 0.1, Es = -0.5))

	putterGrowthEquationBrms  <-  brms::bf(lnMass_g ~ log(exp(lnAsympMass + Es * invKT)) + (1 / (1 - gAlpha)) * log(1 - exp(-exp(lnGoTs + Eg * invKT) * (1 - gAlpha) * timeDays * exp(lnAsympMass + Es * invKT)^(gAlpha - 1))), lnGoTs ~ 1 + (1|G|fishID), lnAsympMass ~ 1 + (1|G|fishID), Eg ~ 1, Es ~ 1, gAlpha ~ 1 + (1|G|fishID), nl = TRUE)

	set.seed(1)
	brms::brm(putterGrowthEquationBrms, data = data, family = gaussian(), prior = putterGrowthPriors, sample_prior = TRUE, chains = 3, cores = 3, iter = 3e4, warmup = 1.5e4, inits = putterGrowthInits, control = list(adapt_delta = 0.99, max_treedepth = 15))
}

boltzmannModelBrms  <-  function (data) {
	set.seed(1)
    boltzInits    <-  list(list(lnMass = 0.4, invKT = 0.8),
                           list(lnMass = 0.8, invKT = 0.75),
                           list(lnMass = 1.0, invKT = 0.6))

	brms::brm(lnRespRate_j_d ~ lnMass_g + invKT + (1 + lnMass_g | fishID), data = data, family = gaussian(), prior = c(prior(normal(1, 2), 'b'), prior(normal(3, 3), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 3, cores = 3, iter = 3e4, warmup = 1.5e4, inits = boltzInits, control = list(adapt_delta = 0.99, max_treedepth = 15))
}

schoolfieldModelBrms  <-  function (data, response) {
	schoolfieldPriors    <-  c(prior(normal(1, 1e6), nlpar = 'lnBoTs'),
	                           prior(normal(1, 1e6), nlpar = 'scalingAlpha'),
	                           prior(normal(2, 2), nlpar = 'Ei', lb = 1e-10, ub = 5),
	                           prior(normal(-1.2, 1), nlpar = 'logitEr'),
	                           prior(normal(298, 5), nlpar = 'Topt', lb = 250, ub = 320))

	schoolfieldInits    <-  list(list(scalingAlpha = 0.8, logitEr = -0.2, lnBoTs = -6, Ei = 1, Topt = 295),
	                             list(scalingAlpha = 0.7, logitEr = -1,   lnBoTs = -7, Ei = 2, Topt = 300),
	                             list(scalingAlpha = 0.6, logitEr = -0.7, lnBoTs = -8, Ei = 3, Topt = 305))

	schoolFieldEquation  <-  brms::bf(lnRespRate_j_d ~ lnBoTs + scalingAlpha * lnMass_g + (Ei / (1 + exp((-1) * logitEr))) * invKT - log(1 + (Ei / (1 + exp((-1) * logitEr))) / (Ei - (Ei / (1 + exp((-1) * logitEr)))) * exp(Ei / boltzmannK * (1 / Topt - 1 / tempKelvin))), lnBoTs ~ 1 + (1|G|fishID), scalingAlpha ~ 1 + (1|G|fishID), logitEr ~ 1, Topt ~ 1, Ei ~ 1, nl = TRUE)

	set.seed(1)
	brms::brm(schoolFieldEquation, data = data, family = gaussian(), prior = schoolfieldPriors, sample_prior = TRUE, chains = 3, cores = 3, iter = 3e4, warmup = 1.5e4, inits = schoolfieldInits, control = list(adapt_delta = 0.99, max_treedepth = 15))
}

poRatioModelBrms  <-  function (data) {
	inits    <-  list(list(mass_g = 3.0, tempCelsius = -0.1),
                      list(mass_g = 1.5, tempCelsius = 0.1),
                      list(mass_g = -3.0, tempCelsius = 0))
	brms::brm(poRatio ~ mass_g + tempCelsius + (1 | tankID), data = data, family = gaussian(), prior = c(prior(normal(1, 2), 'b'), prior(normal(3, 3), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 3, cores = 3, iter = 3e4, warmup = 1.5e4, inits = inits, control = list(adapt_delta = 0.99, max_treedepth = 15))
}

runEmMitoRespModels  <-  function (mitoData) {
	modEmAtp   <-  brms::brm(log(emAtp) ~ invKT, data = mitoData, chains = 3, cores = 3, iter = 1e3)
	modEmOxC   <-  brms::brm(log(emOxC) ~ invKT, data = mitoData[mitoData$tempCelsius > 20, ], chains = 3, cores = 3, iter = 1e3)
	list('modEmAtp' = modEmAtp, 'modEmOxC' = modEmOxC)
}

runAsympMassNo20Model  <-  function (data, model) {
	ranG            <-  coef(model)$fishID
	data            <-  data[!is.na(data$mass_g) & !(data$fishID %in% c('32_2_5', '32_3_2')) & data$tempCelsius != 20, ] # dead on week 2
	data$sE         <-  ranG[data$fishID, 'Estimate', 'Es_Intercept']
	data$asympMass  <-  exp(ranG[data$fishID, 'Estimate', 'lnAsympMass_Intercept'] + data$sE * data$invKT)
	brms::brm(log(asympMass) ~ invKT, data = data, chains = 3, cores = 3, iter = 1e3)
}

docKnitObjectsForRmd  <-  function (cleanedData, respirationData, mitochondrialRespirationData, growthModel, respModelBoltz, poRatioModel, respModelSchool) {
	summaryMass  <-  function (x) {
		data.frame('mean' = LoLinR::rounded(mean(x$mass_g), 2), 'sd' = LoLinR::rounded(sd(x$mass_g), 2))
	}
	cleanBrmsOutput  <-  function (data) {
		rownames(data)  <-  gsub('_Intercept', '', rownames(data))
		colnames(data)  <-  c('estimate', 'se', 'lower', 'upper')
		round(data, 2)
	}
	cleanBrmsFixefOutput  <-  function (model) {
		cleanBrmsOutput(brms::fixef(model))
	}
	cleanBrmsRanefOutput  <-  function (model) {
		cleanBrmsOutput(summary(model)$random[[1]][, 1:4, drop = FALSE])
	}
	massRange           <-  LoLinR::rounded(range(cleanedData$mass_g, na.rm = TRUE), 2)
	massRangeStart      <-  plyr::ddply(cleanedData[cleanedData$Date <= as.Date('2018-05-17'), ], .(tempCelsius), summaryMass)
	massRangeEnd        <-  plyr::ddply(cleanedData[cleanedData$Date > as.Date('2018-07-13'), ], .(tempCelsius), summaryMass)
	growthModelParsFix  <-  cleanBrmsFixefOutput(growthModel)
	respModelParsFix    <-  cleanBrmsFixefOutput(respModelBoltz)
	poModelParsFix      <-  cleanBrmsFixefOutput(poRatioModel)
	poModelParsRan      <-  cleanBrmsRanefOutput(poRatioModel)
	growthModelParsRan  <-  cleanBrmsRanefOutput(growthModel)
	respModelParsRan    <-  cleanBrmsRanefOutput(respModelBoltz)
	schoolModelParsRan  <-  cleanBrmsRanefOutput(respModelSchool)
	schoolErs           <-  calculateErWithCI(respModelSchool)
	em20                <-  exp((respModelParsFix['invKT', 'estimate'] - growthModelParsFix['Eg', 'estimate']) / 8.62e-5 * (1 / 299.65 - 1 / 293.15))
	em32                <-  exp((respModelParsFix['invKT', 'estimate'] - growthModelParsFix['Eg', 'estimate']) / 8.62e-5 * (1 / 299.65 - 1 / 305.15))
	emFold20to32        <-  LoLinR::rounded(em32 / em20, 0)

	list('massRange' = massRange, 'massRangeStart' = massRangeStart, 'massRangeEnd' = massRangeEnd, 'growthModelParsFix' = growthModelParsFix, 'respModelParsFix' = respModelParsFix, 'growthModelParsRan' = growthModelParsRan, 'respModelParsRan' = respModelParsRan, 'emFold20to32' = emFold20to32, 'poModelParsFix' = poModelParsFix, 'schoolModelParsRan' = schoolModelParsRan, 'schoolErs' = schoolErs, 'poModelParsRan' = poModelParsRan)
}

calculateErWithCI  <-  function (brmsModel) {
	posteriors    <-  brms::posterior_samples(brmsModel, params = 'b_')
	posteriorErs  <-  posteriors[, 'b_Ei_Intercept'] / (1 + exp((-1) * posteriors[, 'b_logitEr_Intercept']))
	cis           <-  quantile(posteriorErs, probs = c(0.025, 0.975))
	list('allErs' = posteriorErs, 'mean' = mean(posteriorErs), 'median' = median(posteriorErs), 'lower95CI' = cis[1], 'upper95CI' = cis[2])
}
