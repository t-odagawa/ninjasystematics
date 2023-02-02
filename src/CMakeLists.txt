add_executable(InteractionModelSystematics
	InteractionModelSystematics.cpp)

target_link_libraries(InteractionModelSystematics PRIVATE
	${NINJA_SYST_LIBRARIES}
)

# install the execute in the bin folder
install(TARGETS InteractionModelSystematics DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(InteractionModelMyReWeight
	InteractionModelMyReWeight.cpp)

target_link_libraries(InteractionModelMyReWeight PRIVATE
	${NINJA_SYST_LIBRARIES}
)

# install the execute in the bin folder
install(TARGETS InteractionModelMyReWeight DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(InteractionModelCovariance
	InteractionModelCovariance.cpp)

target_link_libraries(InteractionModelCovariance PRIVATE
	${NINJA_SYST_LIBRARIES}
)

# install the execute in the bin folder
install(TARGETS InteractionModelCovariance DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(DetectorCovariance
	DetectorCovariance.cpp)

target_link_libraries(DetectorCovariance PRIVATE
	${NINJA_SYST_LIBRARIES}
)

# install the execute in the bin folder
install(TARGETS DetectorCovariance DESTINATION "${CMAKE_INSTALL_BINDIR}")