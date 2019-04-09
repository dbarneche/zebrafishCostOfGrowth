tableRounding  <-  function (x, rounding = 2) {
    format(round(x, rounding), nsmall = rounding)
}

makeTableS1  <-  function (model) {
	fixefs  <-  brms::fixef(model)
	ranefs  <-  summary(model)$random$fishID

	pars   <-  c('lnGoTs_Intercept', 'lnAsympMass_Intercept', 'Eg_Intercept', 'Es_Intercept', 'gAlpha_Intercept', 'sd(lnGoTs_Intercept)', 'sd(lnAsympMass_Intercept)', 'sd(gAlpha_Intercept)', 'cor(lnGoTs_Intercept,lnAsympMass_Intercept)', 'cor(lnGoTs_Intercept,gAlpha_Intercept)', 'cor(lnAsympMass_Intercept,gAlpha_Intercept)')

	tableS1  <-  rbind(tableRounding(as.matrix(fixefs[pars[1:5], c('Estimate', 'Q2.5', 'Q97.5')]), 2), 
					   matrix(rep('', 6), 2, 3),
					   tableRounding(as.matrix(ranefs[pars[6:length(pars)], c('Estimate', 'l-95% CI', 'u-95% CI')]), 2)
					  )
	rownames(tableS1)  <-  gsub('_Intercept', '', c(pars[1:5], '', '', pars[6:length(pars)]))
	tableS1
}

makeTableS2  <-  function (model) {
	fixefs  <-  brms::fixef(model)
	ranefs  <-  summary(model)$random$fishID
	
	pars  <-  c('Intercept', 'lnMass_g', 'invKT', 'sd(Intercept)', 'sd(lnMass_g)', 'cor(Intercept,lnMass_g)')

	tableS2  <-  rbind(tableRounding(as.matrix(fixefs[pars[1:3], c('Estimate', 'Q2.5', 'Q97.5')]), 2),
					   matrix(rep('', 6), 2, 3),
					   tableRounding(as.matrix(ranefs[pars[4:6], c('Estimate', 'l-95% CI', 'u-95% CI')]), 2)
					  )
	rownames(tableS2)  <-  c(pars[1:3], '', '', pars[4:6])
	tableS2
}

makeTableS3  <-  function (model) {
	fixefs  <-  brms::fixef(model)
	ranefs  <-  summary(model)$random$fishID

	pars   <-  c('lnBoTs_Intercept', 'scalingAlpha_Intercept', 'logitEr_Intercept', 'Ei_Intercept', 'Topt_Intercept', 'sd(lnBoTs_Intercept)', 'sd(scalingAlpha_Intercept)', 'cor(lnBoTs_Intercept,scalingAlpha_Intercept)')

	tableS3  <-  rbind(tableRounding(as.matrix(fixefs[pars[1:5], c('Estimate', 'Q2.5', 'Q97.5')]), 2), 
					   matrix(rep('', 6), 2, 3),
					   tableRounding(as.matrix(ranefs[pars[6:length(pars)], c('Estimate', 'l-95% CI', 'u-95% CI')]), 2)
					  )
	rownames(tableS3)  <-  gsub('_Intercept', '', c(pars[1:5], '', '', pars[6:length(pars)]))
	tableS3
}

makeTableS4  <-  function (model) {
	fixefs  <-  brms::fixef(model)
	ranefs  <-  summary(model)$random$tankID

	pars   <-  c('Intercept', 'mass_g', 'tempCelsius', 'sd(Intercept)')

	tableS4  <-  rbind(tableRounding(as.matrix(fixefs[pars[1:3], c('Estimate', 'Q2.5', 'Q97.5')]), 2),
					   matrix(rep('', 6), 2, 3),
					   tableRounding(as.matrix(ranefs[pars[4], c('Estimate', 'l-95% CI', 'u-95% CI'), drop = FALSE]), 2)
					  )
	rownames(tableS4)  <-  c(pars[1:3], '', '', pars[4])
	tableS4
}

writeTables  <-  function (dest, table) {
	write.csv(table, dest, row.names = FALSE)
}
