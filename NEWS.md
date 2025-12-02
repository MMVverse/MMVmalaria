# MMVmalaria 1.4.1

# testPKmodelsIQR
* Patched to be compatible with more stinrgent checks in IQRnlmeEst() in IQRtools 2025.05; namely that this function will now remove Tlag1 from the dosing options and modelSpec if test lagtime is set to FALSE. 

# Model library
* 1, 3 and 10 step turnover clearance models (abs0 and abs1) added to model library. 