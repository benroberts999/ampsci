# Writing Modules

\brief Instructions for writing custom ampsci modules

- The modules system allows the easy calculation of any atomic properties after the wavefunction has been calculated.
- Any number of _modules_ can be run by adding a `Module::moduleName{}' block to the input file.
- Get a list of available modules: `$ ./ampsci -m`
- See [doc/modules.md](/doc/modules.md) for details of currently avalable modules
- The code is designed so that you can easily create your own modules.

## Creating your own module

An example module is provided to help you write your own module;

- [src/Modules/exampleModule.hpp](/src/Modules/exampleModule.hpp)
- [src/Modules/exampleModule.cpp](/src/Modules/exampleModule.cpp)

You should duplicate this module (both the .cpp and .hpp files) and give it a new name. That will be much easier than starting from scratch.

- Modules are functions that have the following function signature:

  ```cpp
  namespace Module{
  void exampleModule(const IO::InputBlock &input, const Wavefunction &wf);
  }
  ```

- `input` is an `IO::InputBlock` that holds any input options.

- `wf` is the `Wavefunction` object that was calculated by ampsci

- They are typically placed in the `Module` namespace (but don't need to be).

- Typically, modules live in the `src/Modules` directory, but they can live anywhere

- Then, you can do anything you like inside this function

## Including your Module into ampsci

- In order for ampsci to know about your module, you must update the file `src/Modules/module_list.hpp`
  - Add the corresponding header file to the #includes list

- Add a `std::pair` the the `std::vector` _module_list_ in the form

    ```cpp
    {"moduleName", &moduleName}
    ```

- The first of the pair is a string, which will be the name of the module. This is how you will refer to the module in the input file.

- The second of the pair is the pointer to the function name.

- You then have to recompile ampsci (modules are compiled into ampsci for now)

- That's it - you're now ready to run your module by adding the `Module::moduleName{}` block to the input file.
