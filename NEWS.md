# MMVmalaria 1.2.0

This MMVmalaria release is the final release before migration of MMVverse packages to R 4.4.1 (from R 3.6.3). This release is intended to capture changes made to MMVmalaria (and MMVbase) between the previous releases and the subsequent release post-migration. This version of MMVmalaria has not been extensively tested, and it is adviced to use the previous release (1.1.0) for conducting activities in R 3.6.3. 

# Removed functions 
* adjust_Dose
* get_IQRdistribution
* convert_catCovToTxt
* generate_ActivityContent
* generate_ActivityResults
* generate_DuplicatedActivityReport
* generate_MissingActivityReport
* get_ActivityInfo 
* get_ActivityList
* cbindMMV
* rbindMMV
* combMMV
* libraryMMV
* reverse_List

# Functions migrated to MMVbase 
* MMVggplot
* transform_IQRggplotToMMVggplot
* get_ActivityPath
* adjust_ggplotTheme
* summarizeReplicates
* swapName_MMVnameToName
* swapName_NameToMMVname
* transform_dataFrame_LongtoWide
* transform_dataFrame_WidetoLong
* Check_MissingDatabyNAME
* IQRtableToDataFrame
* aux_CommonSubPath
* aux_createUSUBJID
* aux_removeEscapeChar
* check_dataGeneralMMV
* check_dataGeneralMMV_Center
* clopperPearsonMMV
* collapseMMV
* convert_Unit
* create_PKPDsimtime
* file.copyMMV
* find_MinMMV
* is.fileMMV
* getDoselevel
* get_IXGDFtoRemove
* get_fileNameExtension
* newtonRaphson.Function
* newtonRaphson.Vector
* summarize_PIandCIgeneric
* summaryAcrossTrials
* summaryByTrial
* Load_Rdata

# Functions migrated to MMVworkflows
* check_ActivityInfo 
* check_ActivityResults
* create_RprofileMMV
* generate_OverviewFile
* list_ActivityType
* list_Center
* list_GenericType
* saveActivityInfo
* set_MMVsettings

# DESCRIPTION tidying 
* IQRtools moved to Imports (from Depends)
* MMVbase Imports MMVbase (== 1.1.0) (previously no version) 


 