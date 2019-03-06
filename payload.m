function payload = inputweight()
%Array including names, seat numbers and respective weights

continuee = 1
payload = []
while continuee==1
    ques=input('Add another mob member y/n? : ','s');
    if ques=='y'
        infa=input('Name, seat number (1-10), weight[kg] as column ->');
        payload=[payload;infa]
    else continuee=0
    end
end
end