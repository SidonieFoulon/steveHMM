
source("test_functions.r")
testQNEM <- test.qnem(1:50)
save(testQNEM, file = "testQNEM.rda")

testEM <- test.em(1:50)
save(testEM, file = "testEM.rda")

testSQUAREM <- test.em(1:50)
save(testSQUAREM, file = "testSQUAREM.rda")


