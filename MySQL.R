library("RMySQL")
# Create a connection Object to MySQL database.
# We will connect to the sampel database named "testdb" that comes with MySql installation.
mysqlconnection = dbConnect(MySQL(), user = 'root', password = 'fqyZ+kt(U26Y', dbname = 'prediction-servers',
                            host = 'localhost')

# List the tables available in this database.
dbListTables(mysqlconnection)

a = fetch(dbSendQuery(mysqlconnection, "select * from ensembl_peptide_sequences")
, n = 5)
a

# Use the R data frame "mtcars" to create the table in MySql.
# All the rows of mtcars are taken inot MySql.
dbWriteTable(mysqlconnection, "mtcars", mtcars[, ], overwrite = TRUE)