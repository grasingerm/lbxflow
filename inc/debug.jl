macro NDEBUG()
  return false;
end

macro warn(condition, message)
  if @NDEBUG()
    return;
  else
    return quote
      if condition
        println(message);
      end
    end;
  end
end
