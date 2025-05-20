cmake_minimum_required(VERSION 3.19)

project(WebAssembly LANGUAGES CXX)

message(STATUS "Running cmakelist of project WebAssembly...")

find_package(Qt6 6.5 REQUIRED COMPONENTS Core Widgets)
find_package(fmt CONFIG REQUIRED)
find_package(Boost REQUIRED)

if(Boost_FOUND)
  message(STATUS "Boost_INCLUDE_DIRS = '${Boost_INCLUDE_DIRS}'")
  message(STATUS "Boost_LIBRARIES = '${Boost_LIBRARIES}'")
else()
  message(STATUS "Boost was NOT found")
endif()

qt_standard_project_setup()

qt_add_executable(WebAssembly
    WIN32 MACOSX_BUNDLE
    main.cpp
    mainwindow.cpp
    mainwindow.h
    mainwindow.ui
    vcpkg.json

)

target_link_libraries(WebAssembly
    PRIVATE
        Qt::Core
        Qt::Widgets
)
target_link_libraries(WebAssembly PRIVATE fmt::fmt)
#target_link_libraries(WebAssembly PRIVATE Boost::numeric_odeint)
target_link_libraries(WebAssembly PRIVATE ${Boost_LIBRARIES})


include(GNUInstallDirs)

install(TARGETS WebAssembly
    BUNDLE  DESTINATION .
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

qt_generate_deploy_app_script(
    TARGET WebAssembly
    OUTPUT_SCRIPT deploy_script
    NO_UNSUPPORTED_PLATFORM_ERROR
)
install(SCRIPT ${deploy_script})


