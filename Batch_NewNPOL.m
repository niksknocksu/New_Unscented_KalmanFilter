% Run multiple NPOL scnearios mentioned in NPOL_new.xlsx starting sheet 1
sheet = 1;
MC = 10000;
while sheet <=11
    fprintf('Scenario %d \n', sheet);
    func_NewNPOL_NUKF(sheet,MC)
    sheet = sheet + 1;
end