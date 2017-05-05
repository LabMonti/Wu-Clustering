function ClearAppData(var_name,hObject)
    if isappdata(hObject.Parent,var_name)
        rmappdata(hObject.Parent,var_name)
    end
end
            