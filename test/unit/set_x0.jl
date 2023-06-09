
using BiochemNetABC

load_model("SIR")
load_model("ER")

test_all = true
new_x0 = [92,10,2]
set_x0!(SIR, new_x0)
test_all = test_all && SIR.x0 == new_x0

new_x0 = [1, 1, 4, 2]
set_x0!(ER, new_x0)
test_all = test_all && ER.x0 == new_x0 && SIR.x0 == [92, 10, 2]

new_x0 = [10,10,2]
set_x0!(SIR, new_x0)
test_all = test_all && SIR.x0 == new_x0 && ER.x0 == [1, 1, 4, 2]

return test_all

