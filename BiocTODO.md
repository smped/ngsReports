# DESCRIPTION

- [x] Point your BugReports to the issues page
- [x] Provide a more complete Description field
- [x] Avoid using any GitHub repositories (remove the Remotes field)

# NAMESPACE
- [ ] Stick to either snake_case or camelCase (snake case is recommended for tidyverse interfacing packages, otherwise use camelCase)
- [ ] Consider shortening function names, long function names can be cumbersome and awkward for the user unless they are established file standards
- [ ] Consider using plot methods on classes so that you don't have to append each plot method with plot*

# Vignettes

- [x] Remove the R file from the folder
- [x] Respect the 80 column width limits. It improves readability and maintainership
- [ ] Use system.file to provide the template.Rmd
- [ ] Do you need to make a distinction between a single file and a list of files? Can they be accomodated by a single FastqcFileList class?
- [ ] It's easy to get lost in the 'Generating Plots' section. Perhaps it would be helpful to categorize groups of functions.

# R

- [ ] As mentioned earlier, you can reduce FastqcFile and FastqcFileList into
a single class. This would avoid you the repitition of methods for both
classes, otherwise use inheritance.
- [ ] Group files where possible to reduce the number of files to maintain, such
as helper functions
- [ ] I don't quite understand why you need to create a generic from the class
name instead of a simple constructor function (FastqcFile).
- [ ] Why not use the established path generic rather than creating fileName?
- [ ] Create a coercion method instead of a class method to move from one class
representation to another
- [ ] Only set methods for classes that are your own and use the ANY class for
setting methods for vectors.
- [ ] It looks like you're not taking advantage of existing Bioconductor classes. Your main class in use is the tibble.
- [ ] Is there functionality that you could use from the ShortRead package?