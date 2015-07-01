macro profif(condition, statement)
  quote
    if $condition
      @profile $statement
    else
      $statement
    end
  end
end
