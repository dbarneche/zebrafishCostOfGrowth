######################
# AUXILLIARY FUNCTIONS
######################
toDev  <-  function (expr, dev, filename, ..., verbose = TRUE) {
    if (verbose) {
        cat(sprintf('Creating %s\n', filename))
    }
    dev(filename, ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
}

toPdf  <-  function (expr, filename, ...) {
    toDev(expr, pdf, filename, ...)
}

timeAtAsymp  <-  function (c1, rE, invKT, alpha, asympMass, gE) {
    # time required to reach 95% of asymptotic mass; assume starting mass at birth = 0
    # c1 = bo / Em
    - (1 / (c1 * exp(rE * invKT))) * log(1 - 0.95^(1 - alpha)) / ((asympMass * exp(gE * invKT))^(alpha - 1) * (1 - alpha))
}

massAtTime  <-  function (asympMass, gE, invKT, goTs, rE, alpha, timeDays) {
    asympMass * exp(gE * invKT) * (1 - exp(-goTs * exp(rE * invKT) * (1 - alpha) * timeDays * (asympMass * exp(gE * invKT))^(alpha - 1)))^(1 / (1 - alpha))
}

massOptimum  <-  function (asympMass, alpha) {
    asympMass * (alpha / 1) ^ (1 / (1 - alpha))
}

##########
# ANALYSES
##########
makeFig1  <-  function (dest, ...) {
    toPdf(fig1(...), dest, width = 7, height = 7)
    extrafont::embed_fonts(dest)
}

fig1  <-  function (data, model) {
    data    <-  data[!is.na(data$mass_g) & !(data$fishID %in% c('32_2_5', '32_3_2')), ] # dead on week 2
    tRange  <-  c(25, 350)
    days    <-  seq(min(tRange), max(tRange), length.out = 60)
    out     <-  plyr::daply(data, .(fishID), function (x, model, days) {
        exp(predict(model, newdata = data.frame(invKT = unique(x$invKT), timeDays = days, fishID = unique(x$fishID)), summary = TRUE))[, 'Estimate']
    }, model = model, days = days)
    
    par(omi = rep(0.5, 4), cex = 1, cex.axis = 1.3, cex.lab = 1.4)
    plot(NA, ylab = 'Mass (g)', xlab = 'Days since fertilisation', las = 1, xlim = tRange, ylim = range(out))
    for (i in 1:nrow(out)) {
        lines(days, out[i, ], lty = 1, lwd = 0.5, col = unique(data$color[data$fishID == rownames(out)[i]]))
    }
    labPos  <-  seq(0.95, 0.725, length.out = 4)
    for (k in seq_along(unique(data$invKT))) {
        invKT  <-  sort(unique(data$invKT), decreasing = TRUE)[k]
        LoLinR::proportionalLabel(c(0.03, 0.06), rep(labPos[k], 2), text = FALSE, type = 'l', lty = 1, lwd = 1, col = unique(data$color[data$invKT == invKT]))
        LoLinR::proportionalLabel(0.08, labPos[k], substitute(a*degree*'C', list(a = unique(data$tempCelsius[data$invKT == invKT]))), adj = c(0, 0.5))
    }
    usr  <-  par('usr')
    lines(c(107, 107), c(usr[3], usr[4]), lty = 2, lwd = 0.8, col = 'black')
}

makeFig2  <-  function (dest, ...) {
    toPdf(fig2(...), dest, width = 7, height = 7)
    extrafont::embed_fonts(dest)
}

fig2  <-  function (data, growthModel, asympMassNo20Model) {
    ranG            <-  coef(growthModel)$fishID
    data            <-  data[!is.na(data$mass_g) & !(data$fishID %in% c('32_2_5', '32_3_2')), ] # dead on week 2
    data$sE         <-  ranG[data$fishID, 'Estimate', 'Es_Intercept']
    data$asympMass  <-  exp(ranG[data$fishID, 'Estimate', 'lnAsympMass_Intercept'] + data$sE * data$invKT)

    par(omi = rep(0.5, 4), cex = 1, cex.axis = 1.3, cex.lab = 1.4)
    plot(asympMass ~ invKT, data = data, ylab = '', xlab = substitute('Inverse temperature (eV'^{-1}*')'), las = 1, cex = 1.4, col = LoLinR::transparentColor(data$color, 0.1), bg = LoLinR::transparentColor(data$color, 0.1), pch = data$shape, log = 'y', xlim = 1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + c(19.5, 32.5))), ylim = c(0.17, 1), axes = FALSE)
    axis(1, at = round(1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + seq(20, 32, 4))), 2))
    axis(3, at = round(1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + seq(20, 32, 4))), 2), labels = seq(20, 32, 4))
    axis(2, las = 1)
    box()
    text(data$invKT[data$tempCelsius == 20], mean(data$asympMass[data$tempCelsius == 20]), '*', adj = c(0.5, 0.5), cex = 2)
    LoLinR::proportionalLabel(-0.17, 0.5, substitute('Asymptotic mass, ' * italic('M(T)') * ' (g)'), srt = 90, cex = 1.4, xpd = NA, log = 'y', adj = c(0.5, 0.5))
    LoLinR::proportionalLabel(0.5, 1.17, substitute('Temperature (' * degree * 'C)'), cex = 1.4, xpd = NA, log = 'y', adj = c(0.5, 0.5))

    data2  <-  data[data$tempCelsius > 20, ]
    coefs  <-  brms::fixef(asympMassNo20Model)
    lines(range(data2$invKT), exp(coefs[1, 1] + coefs[2, 1] * range(data2$invKT)), lty = 2, lwd = 2)
    LoLinR::proportionalLabel(0.45, 0.15, substitute('ln '*italic(y) == a + italic(x) %.% b, list(a = round(coefs[1, 1], 2), b = round(coefs[2, 1], 2))), log = 'y', adj = c(0, 0.5), cex = 1.2)
    LoLinR::proportionalLabel(0.45, 0.05, paste0('Slope 95% CI: ', round(coefs[2, 3], 2), ' - ' , round(coefs[2, 4], 2)), log = 'y', adj = c(0, 0.5), cex = 1.2)
}

makeFig3  <-  function (dest, ...) {
    toPdf(fig3(...), dest, width = 8.5, height = 4.5)
    extrafont::embed_fonts(dest)
}

fig3  <-  function (data, model) {
    coefs              <-  coef(model)$fishID
    fixefs             <-  brms::fixef(model)
    data$lnBoTs        <-  coefs[data$fishID, 'Estimate', 'Intercept']
    lnBoTs             <-  fixefs['Intercept', 'Estimate']
    data$scalingAlpha  <-  coefs[data$fishID, 'Estimate', 'lnMass_g']
    scalingAlpha       <-  fixefs['lnMass_g', 'Estimate']
    data$Er            <-  coefs[data$fishID, 'Estimate', 'invKT']
    Er                 <-  fixefs['invKT', 'Estimate']
    data$boltzmannK    <-  8.62e-5
    boltzmannK         <-  8.62e-5

    data$lnRespMassCorrected  <-  data$lnRespRate_j_d - data$scalingAlpha * data$lnMass_g + data$scalingAlpha * mean(data$lnMass_g)
    data$lnRespTempCorrected  <-  data$lnRespRate_j_d - data$Er * data$invKT + data$Er * mean(data$invKT)
    
    par(mfrow = c(1, 2), mai = c(1.2, 1.02, 0.82, 0.40), omi = rep(0, 4), cex = 1, cex.lab = 1.3, mgp = c(3, 0.5, 0), tck = -0.03, cex.axis = 1.1)
    plot(exp(lnRespTempCorrected) ~ exp(lnMass_g), data = data, col = LoLinR::transparentColor(data$color, 0.5), bg = LoLinR::transparentColor(data$color, 0.5), pch = data$shape, xlab = 'Body mass (g)', ylab = '', ylim = c(4, 140), las = 1, cex = 1.2, log = 'xy')
    title(ylab = 'Oxygen consumption rate', cex = 1.3, line = 3.3, xpd = NA)
    title(ylab = substitute('(J / d) @ ' * a * degree * 'C', list(a = round(mean(data$tempCelsius), 1))), cex = 1.3, line = 1.8, xpd = NA)
    LoLinR::proportionalLabel(0.03, 0.92, substitute(italic(y) == a %.% italic(x)^b, list(a = round(exp(round(lnBoTs, 2)), 2), b = round(scalingAlpha, 2))), log = 'xy', adj = c(0, 0.5), cex = 0.9)

    xseq   <-  seq(min(data$lnMass_g), max(data$lnMass_g), length.out = 50)
    preds  <-  predict(model, newdata = data.frame(lnMass_g = xseq, invKT = mean(data$invKT)), re_formula = NA, summary = TRUE, probs = c(0.025, 0.975))
    lines(exp(xseq), exp(preds[, 1]), lty = 2, lwd = 1.5)
    lines(exp(xseq), exp(preds[, 3]), lty = 2, lwd = 0.8, col = 'grey60')
    lines(exp(xseq), exp(preds[, 4]), lty = 2, lwd = 0.8, col = 'grey60')
    labPos  <-  seq(0.3, 0.05, length.out = 4)
    for (k in seq_along(unique(data$invKT))) {
        invKT  <-  sort(unique(data$invKT), decreasing = TRUE)[k]
        LoLinR::proportionalLabel(0.82, labPos[k], text = FALSE, bg = unique(data$color[data$invKT == invKT]), pch = unique(data$shape[data$invKT == invKT]), cex = 1.2, log = 'xy')
        LoLinR::proportionalLabel(0.85, labPos[k], substitute(a*degree*'C', list(a = unique(data$tempCelsius[data$invKT == invKT]))), adj = c(0, 0.5), log = 'xy')
    }

    par(mai = c(1.2, 0.82, 0.82, 0.60))
    plot(exp(lnRespMassCorrected) ~ invKT, data = data, col = LoLinR::transparentColor(data$color, 0.5), bg = LoLinR::transparentColor(data$color, 0.5), pch = data$shape, xlab = substitute('Inverse temperature (eV'^{-1}*')'), ylab = '', cex = 1.2, xlim = 1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + c(19.5, 32.5))), ylim = c(4, 140), log = 'y', axes = FALSE)
    axis(1, at = round(1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + seq(20, 32, 4))), 2))
    axis(3, at = round(1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + seq(20, 32, 4))), 2), labels = seq(20, 32, 4))
    axis(2, las = 1)
    box()
    title(ylab = 'Oxygen consumption rate', cex = 1.3, line = 3.3, xpd = NA)
    title(ylab = substitute('(J / d) @ ' * a * ' g', list(a = round(exp(mean(data$lnMass_g)), 2))), cex = 1.3, line = 1.8, xpd = NA)
    LoLinR::proportionalLabel(0.5, 1.2, substitute('Temperature (' * degree * 'C)'), cex = 1.3, xpd = NA, log = 'y', adj = c(0.5, 0.5))
    LoLinR::proportionalLabel(0.03, 0.92, substitute(italic(y) == a %.% e^{italic(x) %.% b}, list(a = round(exp(round(lnBoTs, 2)), 2), b = round(Er, 2))), log = 'y', adj = c(0, 0.5), cex = 0.9)

    xseq   <-  seq(20, 32, 1)
    preds  <-  predict(model, newdata = data.frame(lnMass_g = mean(data$lnMass_g), invKT = (1 / boltzmannK * (1 / mean(data$tempKelvin) - 1 / (xseq + 273.15)))), re_formula = NA, summary = TRUE, probs = c(0.025, 0.975))
    lines((1 / boltzmannK * (1 / mean(data$tempKelvin) - 1 / (xseq + 273.15))), exp(preds[, 1]), lty = 2, lwd = 1.5)
    lines((1 / boltzmannK * (1 / mean(data$tempKelvin) - 1 / (xseq + 273.15))), exp(preds[, 3]), lty = 2, lwd = 0.8, col = 'grey60')
    lines((1 / boltzmannK * (1 / mean(data$tempKelvin) - 1 / (xseq + 273.15))), exp(preds[, 4]), lty = 2, lwd = 0.8, col = 'grey60')
    
}

makeFig4  <-  function (dest, ...) {
    toPdf(fig4(...), dest, width = 8.5, height = 4.5)
    extrafont::embed_fonts(dest)
}

fig4  <-  function (data, model) {
    coefs            <-  coef(model)$tankID
    fixefs           <-  brms::fixef(model)
    data$intercept   <-  coefs[data$tankID, 'Estimate', 'Intercept']
    intercept        <-  fixefs['Intercept', 'Estimate']
    data$massSlope   <-  coefs[data$tankID, 'Estimate', 'mass_g']
    massSlope        <-  fixefs['mass_g', 'Estimate']
    data$tempSlope   <-  coefs[data$tankID, 'Estimate', 'tempCelsius']
    tempSlope        <-  fixefs['tempCelsius', 'Estimate']

    data$respMassCorrected  <-  data$poRatio - data$massSlope * data$mass_g + data$massSlope * mean(data$mass_g)
    data$respTempCorrected  <-  data$poRatio - data$tempSlope * data$tempCelsius + data$tempSlope * mean(data$tempCelsius)
    
    par(mfrow = c(1, 2), mai = c(1.2, 1.02, 0.82, 0.40), omi = rep(0, 4), cex = 1, cex.lab = 1.3, mgp = c(3, 0.5, 0), tck = -0.03, cex.axis = 1.1)
    plot(respTempCorrected ~ mass_g, data = data, col = LoLinR::transparentColor(data$color, 0.5), bg = LoLinR::transparentColor(data$color, 0.5), pch = data$shape, xlab = 'Body mass (g)', ylab = '', xlim = c(0.1, 0.6), ylim = c(0, 2), las = 1, cex = 1.2)
    title(ylab = substitute('P/O ratio @ ' * a * degree * 'C', list(a = round(mean(data$tempCelsius), 1))), cex = 1.3, line = 2.5, xpd = NA)
    LoLinR::proportionalLabel(0.03, 0.92, substitute(italic(y) == a + b %.% italic(x), list(a = round(intercept, 2), b = round(massSlope, 2))), adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.03, 0.84, paste0('Slope 95% CI: ', round(fixefs['mass_g', 3], 2), ' - ' , round(fixefs['mass_g', 4], 2)), adj = c(0, 0.5), cex = 0.9)

    xseq   <-  seq(min(data$mass_g), max(data$mass_g), length.out = 50)
    preds  <-  predict(model, newdata = data.frame(mass_g = xseq, tempCelsius = mean(data$tempCelsius)), re_formula = NA, summary = TRUE, probs = c(0.025, 0.975))
    lines(xseq, preds[, 1], lty = 2, lwd = 1.5)
    lines(xseq, preds[, 3], lty = 2, lwd = 0.8, col = 'grey60')
    lines(xseq, preds[, 4], lty = 2, lwd = 0.8, col = 'grey60')

    par(mai = c(1.2, 0.82, 0.82, 0.60))
    plot(respMassCorrected ~ tempCelsius, data = data, col = LoLinR::transparentColor(data$color, 0.5), bg = LoLinR::transparentColor(data$color, 0.5), pch = data$shape, xlab = substitute('Temperature ('*degree*'C)'), ylab = '', cex = 1.2, xlim = c(19.5, 32.5), ylim = c(0, 2), las = 1)
    title(ylab = substitute('P/O ratio @ ' * a * ' g', list(a = round(mean(data$mass_g), 2))), cex = 1.3, line = 2.5, xpd = NA)
    LoLinR::proportionalLabel(0.03, 0.92, substitute(italic(y) == a + b %.% italic(x), list(a = round(intercept, 2), b = round(tempSlope, 2))), adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.03, 0.84, paste0('Slope 95% CI: ', round(fixefs['tempCelsius', 3], 2), ' - ' , round(fixefs['tempCelsius', 4], 2)), adj = c(0, 0.5), cex = 0.9)

    xseq   <-  seq(20, 32, 1)
    preds  <-  predict(model, newdata = data.frame(mass_g = mean(data$mass_g), tempCelsius = xseq), re_formula = NA, summary = TRUE, probs = c(0.025, 0.975))
    lines(xseq, preds[, 1], lty = 2, lwd = 1.5)
    lines(xseq, preds[, 3], lty = 2, lwd = 0.8, col = 'grey60')
    lines(xseq, preds[, 4], lty = 2, lwd = 0.8, col = 'grey60')
}

makeFig5  <-  function (dest, ...) {
    toPdf(fig5(...), dest, width = 7, height = 7)
    extrafont::embed_fonts(dest)
}

fig5  <-  function (data, respModel, growthModel) {
    data            <-  data[!is.na(data$mass_g) & !(data$fishID %in% c('32_2_5', '32_3_2')), ] # dead on week 2
    ranR            <-  coef(respModel)$fishID
    ranG            <-  coef(growthModel)$fishID
    hist(ranR[unique(data$fishID), 'Estimate', 'lnMass_g'] - ranG[unique(data$fishID), 'Estimate', 'gAlpha_Intercept'], col = 'grey', xlab = substitute(italic(alpha) - italic(theta)), las = 1, ylim = c(0, 22), main = '')
}

makeFig6  <-  function (dest, ...) {
    toPdf(fig6(...), dest, width = 8.5, height = 4.5)
    extrafont::embed_fonts(dest)
}

fig6  <-  function (data, emMitoRespModels) {
    par(mfrow = c(1, 2), mai = c(1.2, 0.92, 0.82, 0.50), omi = rep(0, 4), cex = 1, cex.lab = 1.3, mgp = c(3, 0.5, 0), tck = -0.03, cex.axis = 1.1)
    plot(emAtp ~ invKT, data = data, col = LoLinR::transparentColor(data$color, 0.5), bg = LoLinR::transparentColor(data$color, 0.5), pch = data$shape, xlab = '', ylab = substitute('Cost of growth,'~italic(E[m])*', (mM ATP / g)'), ylim = c(0.5, 50), xlim = 1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + c(19.5, 32.5))), log = 'y', cex = 1.2, axes = FALSE)
    axis(1, at = round(1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + seq(20, 32, 4))), 2))
    axis(3, at = round(1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + seq(20, 32, 4))), 2), labels = seq(20, 32, 4))
    axis(2, las = 1)
    box()
    LoLinR::proportionalLabel(0.5, 1.19, substitute('Temperature (' * degree * 'C)'), cex = 1.3, xpd = NA, log = 'y', adj = c(0.5, 0.5))
    LoLinR::proportionalLabel(0.5, -0.19, substitute('Inverse temperature (eV'^{-1}*')'), cex = 1.3, xpd = NA, log = 'y', adj = c(0.5, 0.5))
    coefs  <-  brms::fixef(emMitoRespModels$modEmAtp)
    lines(range(data$invKT), exp(coefs[1, 1] + coefs[2, 1] * range(data$invKT)), lty = 2, lwd = 2)
    LoLinR::proportionalLabel(0.03, 0.92, substitute('ln '*italic(y) == a + italic(x) %.% b, list(a = round(coefs[1, 1], 2), b = round(coefs[2, 1], 2))), log = 'y', adj = c(0, 0.5), cex = 1)
    LoLinR::proportionalLabel(0.03, 0.84, paste0('Slope 95% CI: ', round(coefs[2, 3], 2), ' - ' , round(coefs[2, 4], 2)), log = 'y', adj = c(0, 0.5), cex = 1)

    par(mai = c(1.2, 0.62, 0.82, 0.80))
    plot(emOxC ~ invKT, data = data, col = LoLinR::transparentColor(data$color, 0.5), bg = LoLinR::transparentColor(data$color, 0.5), pch = data$shape, ylab = '', xlab = '', ylim = c(1, 10), xlim = 1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + c(19.5, 32.5))), log = 'y', cex = 1.2, axes = FALSE, xpd = NA)
    title(ylab = substitute('Cost of growth,'~italic(E[m])*', (kJ / g)'), cex.lab = 1.3, line = 2, xpd = NA)
    LoLinR::proportionalLabel(0.5, 1.19, substitute('Temperature (' * degree * 'C)'), cex = 1.3, xpd = NA, log = 'y', adj = c(0.5, 0.5))
    LoLinR::proportionalLabel(0.5, -0.19, substitute('Inverse temperature (eV'^{-1}*')'), cex = 1.3, xpd = NA, log = 'y', adj = c(0.5, 0.5))
    axis(1, at = round(1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + seq(20, 32, 4))), 2))
    axis(3, at = round(1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + seq(20, 32, 4))), 2), labels = seq(20, 32, 4))
    axis(2, las = 1)
    box()
    coefs  <-  brms::fixef(emMitoRespModels$modEmOxC)
    lines(range(data$invKT[data$tempCelsius > 20]), exp(coefs[1, 1] + coefs[2, 1] * range(data$invKT[data$tempCelsius > 20])), lty = 2, lwd = 2)
    LoLinR::proportionalLabel(0.03, 0.92, substitute('ln '*italic(y) == a + italic(x) %.% b, list(a = round(coefs[1, 1], 2), b = LoLinR::rounded(coefs[2, 1], 2))), log = 'y', adj = c(0, 0.5), cex = 1)
    LoLinR::proportionalLabel(0.03, 0.84, paste0('Slope 95% CI: ', round(coefs[2, 3], 2), ' - ' , round(coefs[2, 4], 2)), log = 'y', adj = c(0, 0.5), cex = 1)
}

timeAtAsymp  <-  function (aoTs, theta, asympMass) {
    # time required to reach 95% of asymptotic mass; assume starting mass at birth = 0
    # aoTs = bo / Em
    - (1 / aoTs) * log(1 - 0.95 ^ (1 - theta)) / (asympMass ^ (theta - 1) * (1 - theta))
}

massAtAgePlot  <-  function (dat, xaxt, yaxt, xlim, preds, days) {
    plot(mass_g ~ timeDays, data = dat, ylab = '', xlab = '', las = 1, type = 'n', xlim = xlim, ylim = c(0, 1), axes = FALSE)
    if (xaxt == 's') {
        axis(1, mgp = c(3, 0.1, 0), tck = -0.03)
    } else {
        axis(1, tck = -0.01, labels = NA)
    }
    if (yaxt == 's') {
        axis(2, mgp = c(3, 0.25, 0), tck = -0.03, las = 1, at = seq(0.1, 0.9, 0.4))
    } else {
        axis(2, tck = -0.03, at = seq(0.1, 0.9, 0.4), labels = NA)
    }
    box()
    polygon(c(days[days <= max(dat$timeDays)], rev(days[days <= max(dat$timeDays)])), c(preds[days <= max(dat$timeDays), 'Q2.5'], rev(preds[days <= max(dat$timeDays), 'Q97.5'])), border = NA, col = LoLinR::transparentColor('grey60', 0.5))
    polygon(c(days[days > max(dat$timeDays)], rev(days[days > max(dat$timeDays)])), c(preds[days > max(dat$timeDays), 'Q2.5'], rev(preds[days > max(dat$timeDays), 'Q97.5'])), border = NA, col = LoLinR::transparentColor('grey30', 0.5))
    points(mass_g ~ timeDays, data = dat, cex = 0.6, bg = dat$color, pch = dat$shape)
    lines(days[days <= max(dat$timeDays)], preds[days <= max(dat$timeDays), 1], lty = 2)
    lines(days[days > max(dat$timeDays)], preds[days > max(dat$timeDays), 1], lty = 2, col = unique(dat$color), lwd = 2)
    LoLinR::proportionalLabel(0.03, 0.95, substitute(italic('T') == a * degree * 'C', list(a = unique(dat$tempCelsius))), cex = 0.6)
    LoLinR::proportionalLabel(0.03, 0.85, paste0('Tank # = ', unique(dat$tankN)), cex = 0.6)
    LoLinR::proportionalLabel(0.03, 0.75, paste0('Fish # = ', unique(dat$Basket)), cex = 0.6)
    LoLinR::proportionalLabel(0.03, 0.65, substitute(italic('M') == a * ' g', list(a = LoLinR::rounded(unique(dat$asympMass), 2))), cex = 0.6)
}

makeFigS1  <-  function (dest, ...) {
    toPdf(figS1(...), dest, width = 8, height = 10)
    extrafont::embed_fonts(dest)
}
 
figS1  <-  function (data, growthModel, temperatures) {
    data               <-  data[!is.na(data$mass_g) & !(data$fishID %in% c('32_2_5', '32_3_2')) & data$tempCelsius %in% temperatures, ] # dead on week 2
    brmsPars        <-  coef(growthModel)$fishID
    data$theta      <-  brmsPars[data$fishID, 'Estimate', 'gAlpha_Intercept']
    data$Es         <-  brmsPars[data$fishID, 'Estimate', 'Es_Intercept']
    data$Eg         <-  brmsPars[data$fishID, 'Estimate', 'Eg_Intercept']
    data$asympMass  <-  exp(brmsPars[data$fishID, 'Estimate', 'lnAsympMass_Intercept'] + data$Es * data$invKT)
    data$aoTs       <-  exp(brmsPars[data$fishID, 'Estimate', 'lnGoTs_Intercept'] + data$Eg * data$invKT)

    par(mfcol = c(8, 5), mai = c(0.0132, 0.1452, 0.0132, 0.0132), omi = c(0.528, 0.628, 0.264, 0.1), cex = 1, cex.axis = 0.7, cex.lab = 1.3)
    maxAsymp  <-  exp(brmsPars[data$fishID, 'Q97.5', 'lnAsympMass_Intercept'] + data$Es * data$invKT)
    xlim      <-  c(min(data$timeDays), timeAtAsymp(aoTs = data$aoTs[which.max(maxAsymp)], theta = data$theta[which.max(maxAsymp)], asympMass = maxAsymp[which.max(maxAsymp)]))

    days     <-  seq(min(xlim), max(xlim), length.out = 60)
    newData  <-  plyr::ddply(data, .(fishID), function (z, days) {
        data.frame(invKT = unique(z$invKT), timeDays = days)
    }, days = days)
    preds     <-  exp(predict(growthModel, newdata = newData, summary = TRUE, probs = c(0.025, 0.975)))
    allFish   <-  sort(unique(data$fishID))
    n         <-  0 #start counting
    for (i in seq_along(allFish)) {
        n  <-  n + 1
        yaxt  <-  ifelse(n %in% 1:8, 's', 'n')
        xaxt  <-  ifelse(n %in% seq(8, 40, 8), 's', 'n')
        x     <-  data[data$fishID == allFish[i], ]
        y     <-  preds[newData$fishID == allFish[i], ]
        massAtAgePlot(dat = x, xaxt = xaxt, yaxt = yaxt, xlim = xlim, preds = y, days = days)
    }
    mtext('Time since fertilisation', side = 1, line = 1.5, outer = TRUE, cex = 1.3)
    mtext('Mass (g)', side = 2, line = 1.2, outer = TRUE, cex = 1.3)
}

makeFigS2  <-  function (dest, ...) {
    toPdf(figS2(...), dest, width = 8, height = 10)
    extrafont::embed_fonts(dest)
}
 
figS2  <-  function (data, respModelBoltz, temperatures) {
    brmsPars    <-  coef(respModelBoltz)$fishID
    data        <-  data[data$tempCelsius %in% temperatures, ]
    data$alpha  <-  brmsPars[data$fishID, 'Estimate', 'lnMass_g']
    data$Er     <-  brmsPars[data$fishID, 'Estimate', 'invKT']
    data$boTs   <-  exp(brmsPars[data$fishID, 'Estimate', 'Intercept'] + data$Er * data$invKT)

    par(mfcol = c(8, 5), mai = c(0.0132, 0.1452, 0.0132, 0.0132), omi = c(0.528, 0.628, 0.264, 0.1), cex = 1, cex.axis = 0.7, cex.lab = 1.3)

    newData  <-  plyr::ddply(data, .(fishID), function (z) {
        data.frame(invKT = unique(z$invKT), lnMass_g = range(z$lnMass_g))
    })
    preds     <-  cbind(newData, predict(respModelBoltz, newdata = newData, summary = TRUE, probs = c(0.025, 0.975)))
    allFish   <-  sort(unique(data$fishID))
    n         <-  0 # start counting
    for (i in seq_along(allFish)) {
        n     <-  n + 1
        yaxt  <-  ifelse(n %in% 1:8, 's', 'n')
        xaxt  <-  ifelse(n %in% seq(8, 40, 8), 's', 'n')
        x     <-  data[data$fishID == allFish[i], ]
        y     <-  preds[preds$fishID == allFish[i], ]
        metRatePlot(dat = x, xaxt = xaxt, yaxt = yaxt, preds = y)
    }
    mtext('Body mass (g)', side = 1, line = 1.5, outer = TRUE, cex = 1.3)
    mtext('Oxygen consumption rate (J / d)', side = 2, line = 1.3, outer = TRUE, cex = 1.3)
}

metRatePlot  <-  function (dat, xaxt, yaxt, preds) {
    plot(exp(lnRespRate_j_d) ~ exp(lnMass_g), data = dat, ylab = '', xlab = '', las = 1, type = 'n', xlim = c(0.03, 0.5), ylim = c(4, 240), axes = FALSE, log = 'xy')
    if (xaxt == 's') {
        axis(1, mgp = c(3, 0.1, 0), tck = -0.03)
    } else {
        axis(1, tck = -0.01, labels = NA)
    }
    if (yaxt == 's') {
        axis(2, mgp = c(3, 0.25, 0), tck = -0.03, las = 1, at = c(5, 10, 20, 50, 100, 200))
    } else {
        axis(2, tck = -0.03, at = c(5, 10, 20, 50, 100, 200), labels = NA)
    }
    box()
    polygon(exp(c(preds$lnMass_g, rev(preds$lnMass_g))), exp(c(preds$Q2.5, rev(preds$Q97.5))), border = NA, col = LoLinR::transparentColor('grey60', 0.5))
    points(exp(lnRespRate_j_d) ~ exp(lnMass_g), data = dat, cex = 0.6, bg = dat$color, pch = dat$shape)
    lines(exp(preds$lnMass_g), exp(preds$Estimate), lty = 2)
    LoLinR::proportionalLabel(0.03, 0.95, substitute(italic('T') == a * degree * 'C', list(a = unique(dat$tempCelsius))), cex = 0.6, log = 'xy')
    LoLinR::proportionalLabel(0.03, 0.85, paste0('Tank # = ', unique(dat$tankN)), cex = 0.6, log = 'xy')
    LoLinR::proportionalLabel(0.03, 0.75, paste0('Fish # = ', unique(dat$Basket)), cex = 0.6, log = 'xy')
    LoLinR::proportionalLabel(0.95, 0.1, substitute(italic(y) == a %.% italic(x)^b, list(a = LoLinR::rounded(unique(dat$boTs), 2), b = LoLinR::rounded(unique(dat$alpha), 2))), log = 'xy', adj = c(1, 0.5), cex = 0.7)
}

makeFigS3  <-  function (dest, ...) {
    toPdf(figS3(...), dest, width = 8.5, height = 4.5)
    extrafont::embed_fonts(dest)
}

figS3  <-  function (data, model) {
    coefs              <-  coef(model)$fishID
    fixefs             <-  brms::fixef(model)
    data$lnBoTs        <-  coefs[data$fishID, 'Estimate', 'lnBoTs_Intercept']
    lnBoTs             <-  fixefs['lnBoTs_Intercept', 'Estimate']
    data$scalingAlpha  <-  coefs[data$fishID, 'Estimate', 'scalingAlpha_Intercept']
    scalingAlpha       <-  fixefs['scalingAlpha_Intercept', 'Estimate']
    data$Ei            <-  coefs[data$fishID, 'Estimate', 'Ei_Intercept']
    Ei                 <-  fixefs['Ei_Intercept', 'Estimate']
    data$logitEr       <-  coefs[data$fishID, 'Estimate', 'logitEr_Intercept']
    logitEr            <-  fixefs['logitEr_Intercept', 'Estimate']
    data$Er            <-  data$Ei / (1 + exp(-1 * data$logitEr))
    Er                 <-  Ei / (1 + exp(-1 * logitEr))
    data$Topt          <-  coefs[data$fishID, 'Estimate', 'Topt_Intercept']
    Topt               <-  fixefs['Topt_Intercept', 'Estimate']
    data$boltzmannK    <-  8.62e-5
    boltzmannK         <-  8.62e-5

    data$lnRespMassCorrected  <-  data$lnRespRate_j_d - data$scalingAlpha * data$lnMass_g + data$scalingAlpha * mean(data$lnMass_g)
    data$lnRespTempCorrected  <-  data$lnRespRate_j_d - data$Er / data$boltzmannK * (1 / mean(data$tempKelvin) - 1 / data$tempKelvin) + log(1 + data$Er / (data$Ei - data$Er) * exp(data$Ei / data$boltzmannK * (1 / data$Topt - 1 / data$tempKelvin)))
    
    par(mfrow = c(1, 2), mai = c(1.2, 1.02, 0.82, 0.40), omi = rep(0, 4), cex = 1, cex.lab = 1.3, mgp = c(3, 0.5, 0), tck = -0.03, cex.axis = 1.1)
    plot(exp(lnRespTempCorrected) ~ exp(lnMass_g), data = data, col = LoLinR::transparentColor(data$color, 0.5), bg = LoLinR::transparentColor(data$color, 0.5), pch = data$shape, xlab = 'Body mass (g)', ylab = '', ylim = c(4, 140), las = 1, cex = 1.2, log = 'xy')
    title(ylab = 'Oxygen consumption rate', cex = 1.3, line = 3.3, xpd = NA)
    title(ylab = substitute('(J / d) @ ' * a * degree * 'C', list(a = round(mean(data$tempCelsius), 1))), cex = 1.3, line = 1.8, xpd = NA)
    LoLinR::proportionalLabel(0.03, 0.92, substitute(italic(y) == a %.% italic(x)^b, list(a = round(exp(round(lnBoTs, 2)), 2), b = round(scalingAlpha, 2))), log = 'xy', adj = c(0, 0.5), cex = 0.9)

    xseq   <-  seq(min(data$lnMass_g), max(data$lnMass_g), length.out = 50)
    preds  <-  predict(model, newdata = data.frame(lnMass_g = xseq, invKT = mean(data$invKT), tempKelvin = mean(data$tempKelvin), boltzmannK = boltzmannK), re_formula = NA, summary = TRUE, probs = c(0.025, 0.975))
    lines(exp(xseq), exp(preds[, 1]), lty = 2, lwd = 1.5)
    lines(exp(xseq), exp(preds[, 3]), lty = 2, lwd = 0.8, col = 'grey60')
    lines(exp(xseq), exp(preds[, 4]), lty = 2, lwd = 0.8, col = 'grey60')

    labPos  <-  seq(0.3, 0.05, length.out = 4)
    for (k in seq_along(unique(data$invKT))) {
        invKT  <-  sort(unique(data$invKT), decreasing = TRUE)[k]
        LoLinR::proportionalLabel(0.82, labPos[k], text = FALSE, bg = unique(data$color[data$invKT == invKT]), pch = unique(data$shape[data$invKT == invKT]), cex = 1.2, log = 'xy')
        LoLinR::proportionalLabel(0.85, labPos[k], substitute(a*degree*'C', list(a = unique(data$tempCelsius[data$invKT == invKT]))), adj = c(0, 0.5), log = 'xy')
    }

    par(mai = c(1.2, 0.82, 0.82, 0.60))
    plot(exp(lnRespMassCorrected) ~ invKT, data = data, col = LoLinR::transparentColor(data$color, 0.5), bg = LoLinR::transparentColor(data$color, 0.5), pch = data$shape, xlab = substitute('Inverse temperature (eV'^{-1}*')'), ylab = '', cex = 1.2, xlim = 1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + c(19.5, 32.5))), ylim = c(4, 140), log = 'y', axes = FALSE)
    axis(1, at = round(1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + seq(20, 32, 4))), 2))
    axis(3, at = round(1 / 8.62e-5 * (1 / mean(data$tempKelvin) - 1 / (273.15 + seq(20, 32, 4))), 2), labels = seq(20, 32, 4))
    axis(2, las = 1)
    box()
    LoLinR::proportionalLabel(0.03, 0.25, substitute(italic('E'['r']) == a * ' eV', list(a = round(Er, 2))), log = 'y', adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.03, 0.18, substitute(italic('E'['i']) == a * ' eV', list(a = round(Ei, 2))), log = 'y', adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.03, 0.11, substitute(italic('T'['opt']) == a * ' ' * degree * 'C', list(a = round(Topt - 273.15, 2))), log = 'y', adj = c(0, 0.5), cex = 0.9)

    title(ylab = 'Oxygen consumption rate', cex = 1.3, line = 3.3, xpd = NA)
    title(ylab = substitute('(J / d) @ ' * a * ' g', list(a = round(exp(mean(data$lnMass_g)), 2))), cex = 1.3, line = 1.8, xpd = NA)
    LoLinR::proportionalLabel(0.5, 1.2, substitute('Temperature (' * degree * 'C)'), cex = 1.3, xpd = NA, log = 'y', adj = c(0.5, 0.5))

    xseq   <-  seq(20, 32, 1)
    preds  <-  predict(model, newdata = data.frame(lnMass_g = mean(data$lnMass_g), invKT = (1 / boltzmannK * (1 / mean((xseq + 273.15)) - 1 / (xseq + 273.15))), tempKelvin = xseq + 273.15, boltzmannK = boltzmannK), re_formula = NA, summary = TRUE, probs = c(0.025, 0.975))
    lines((1 / boltzmannK * (1 / mean(data$tempKelvin) - 1 / (xseq + 273.15))), exp(preds[, 1]), lty = 2, lwd = 1.5)
    lines((1 / boltzmannK * (1 / mean(data$tempKelvin) - 1 / (xseq + 273.15))), exp(preds[, 3]), lty = 2, lwd = 0.8, col = 'grey60')
    lines((1 / boltzmannK * (1 / mean(data$tempKelvin) - 1 / (xseq + 273.15))), exp(preds[, 4]), lty = 2, lwd = 0.8, col = 'grey60')
}
