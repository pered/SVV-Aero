function payloadadd = inputweight()
%Array including names, seat numbers and respective weights

continuee = 1
payloadadd = []
while continuee==1
    ques=input('Add another mobster y/n? : ','s');
    if ques=='y'
        infa=input('Seat number (1-10), weight [lbs] as row ->');
        payloadadd=[payloadadd;infa]
    else continuee=0
    end
end
return